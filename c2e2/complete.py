
import sympy
import numpy as np
import os
import json
import shlex, subprocess
import time

from collections import defaultdict

from c2e2.cpp_backend.lib.libc2e2 import IntegerVector, DoubleVector, Model

from .SymEq import SymEq
from .automaton import Automaton, Mode, Invariant, Transition
from .jacobiancalc import JacobianCalc as jcalc
from .property import Property
from .set import RectangleSet
from .constants import *

class Complete_Automata:
    '''
    TODO: Add definition for simulator
    '''

    # path to the current file(simulate.py) 
    cur_path = os.path.dirname(os.path.abspath(__file__)) + '/'

    def __init__(self, automata=None, initial_set=None, unsafe_set=None, time_step=0.01,time_bound=0, prop=None, name="default simulation", simulator=ODEINT_FIX):
        # force automata to be a list
        if isinstance(automata, Automaton):
            automata = [automata,]
        # NOTE: it's list of automaton
        self.automata = automata
        self.initial_set = initial_set
        self.unsafe_set = unsafe_set
        self.property = prop

        # name for this current simulation
        self.name = name

        # time
        self.time_step = time_step
        self.time_horizon = time_bound

        # refinement strategy default
        self.refine_strat = DEF_STRAT
        # bloated strategy for linear mode
        self.linear_strat = GSTAR

        # Simulator Type {ODEINT_FIX | ODEINT_ADP | CAPD}
        self.simulator = simulator

        # populate Inv Guards
        self.guardResets = None
        self.invariants = None
        
        # simple parse
        self.simple_parse()

    def __repr__(self):
        DISPLAY_INFO = "=== Simulator() ===\n"
        DISPLAY_INFO += "Automata: {}\n".format(self.automata_names)
        return DISPLAY_INFO

    def __call__(self, automata, initial_set=None, unsafe_set=None, time_step=0.01,time_bound=0, prop=None):
        """ update/assign automaton by calling Simulator object directly """
        if isinstance(automata, Automaton):
            automata = [automata,]
        # should be list of automaton
        self.automata = automata
        self.property = prop
        self.initial_set = initial_set
        self.unsafe_set = unsafe_set
        # simple parse
        self.simple_parse()


    def simple_parse(self):
        # check for property
        if self.property is not None:
            if isinstance(self.property, Property):
                self.initial_set = self.property.initial_set
                self.unsafe_set = self.property.unsafe_set
            else:
                raise Exception("[Simulator Property Error]: Input property must be type Property()")
        # check for initial set
        if self.initial_set is not None and not isinstance(self.initial_set, RectangleSet):
            raise Exception("[Simulator Initial Set Error]: Input Initial Set must be type RectangleSet()")
        # check for unsafe set
        if self.unsafe_set is not None and not isinstance(self.unsafe_set, RectangleSet):
            raise Exception("[Simulator Unsafe Set Error]: Input Unsafe Set must be type RectangleSet()")

        # TODO: If both property & set are given, check if they are equal(or just use input initial/unsafe set)
 
    def simulate(self):
        # Compose all
        self.compose_all()

        # Generate simulator.cpp at run-time
        path = '../work-dir/simulator.cpp'
        step_type = 'constant'  # or 'adaptive'
        self.generate_simulator(path, self.automata, step_type=step_type)

        # Generate Invariant & Guard Checking
        self.printHybridSimGuardsInvariants()
        self.printBloatedSimGuardsInvariants()

        #TODO * Load Config *
        with open(Simulator.cur_path + '../config.json') as f:
            config = json.load(f)
        CPP_MODEL = self.initialize_cpp_model(config, SIMULATE)
        self.compile_executable(config)

        # Simulate selected model
        print("Running simulate...\n")
        start_time = time.time()
        result = CPP_MODEL.simulate_verify()
        print("--- " + str(time.time() - start_time) + " seconds ---" )
        
        result_str_dict = {
            1 : "Safe",
            0 : "Unknown",
            -1 : "Unsafe"
        }
        print("RESULT: " + result_str_dict[result] + "\n")


        print('Simulation Finished!')
        return

# ==========================================================================================
# ====== TODO: Mostly from C2E2 frontend, need to check more into the code detail ==========
# ==========================================================================================

    @staticmethod
    def generate_simulator(file_path, automata, **kwargs):
        """
        Generate simulator.cpp at run-time.

        file_path   - the path where simulator.cpp is generated
        automata    - list of Automaton() (or can be single Automaton())
        **kwargs    - step_type and other parameters
        """
        # Get kwargs
        step_type = kwargs.pop('step_type', 'adaptive')

        # Setting spefic variables
        if step_type == 'constant':
            integrator = ['runge_kutta4<state_t> stepper;',
                        'size_t steps = integrate_const(stepper, rhs[cur_mode], '
                        'x, ts, te, dt,'
                        '\tIntObs(trace, times));']
        else:
            integrator = ['auto stepper = make_controlled(abs_err, rel_err, '
                        'runge_kutta_dopri5<state_t>());',
                        'size_t steps = integrate_adaptive(stepper, '
                        'rhs[cur_mode], x, ts, te, dt,',
                        '\tIntObs(trace, times));']

        # if user input single automaton, make it into a list
        if isinstance(automata, Automaton):
            automata = [automata,]

        # Obtain and parse variables
        vars = []
        # for each automaton
        for automaton in automata:
            for var in automaton.vars:
                if (var.scope == 'LOCAL_DATA' and 'dot' not in var.name and 
                        'clock' not in var.name):
                    vars.append(var.name)


        # Iterate through all the mode of the automata to
        # obtain and parse differential algebraic equations
        # dxdt is a list modes, each with a list of equations
        modes = []
        dxdt = []
        del_list = []

        for automaton in automata:
            for i, cur_mode in enumerate(automaton.modes):
                modes.append(cur_mode.name) 
                dxdt.append([])
                orig_eqns = []

                # Find variables with '_dot' and extract rhs
                for dai in cur_mode.dais:
                    if '_dot' in dai.raw:
                        # Split the equation and get lhs index
                        # lhs, rhs = dai.raw.split('=')
                        lhs, rhs = str(dai.expr.lhs), dai.expr.rhs
                        lhs = lhs.split('_dot')[0] 
                        lhs_idx = vars.index(lhs)
                        #print "orig: ", rhs
                        #rhs = SymEq.convert_pow(rhs)
                        #print "convert: ", rhs
                        rhs = str(rhs)
                        # Generate jacobian in correct order
                        orig_eqns.insert(lhs_idx,rhs)

                        # Replace variables with 'x[i]'
                        rhs = dai.expr.rhs
                        for j, var in enumerate(vars):
                            rhs = rhs.subs(var, sympy.Symbol('x['+str(j)+']'))
                            # rhs = re.sub(r'\b%s\b' % var, 'x[' + str(j) + ']', rhs)
                        rhs = SymEq.convert_pow(rhs)


                        # Generate dxdt in correct order
                        dxdt[i].insert(lhs_idx, rhs) 

                # Bloat factor calculation/generation
                var_str = ','.join(vars)
                eqn_str = ','.join(orig_eqns)
                del_elem = jcalc.jacobian(var_str, orig_eqns, i)
                del_list.append(del_elem)

        # Generate python file for bloating
        jcalc.createCDFfunction(del_list)

        # Generate dxdt
        for i, mode in enumerate(modes):
            for j, eqn in enumerate(dxdt[i]):
                dxdt[i][j] = 'dxdt[' + str(j) + ']=' + eqn
            dxdt[i] = '\t\t\t' + ';\n\t\t\t'.join(dxdt[i]) + ';'

        # Components for .cpp
        # Auto generated comment
        auto_gen = '/* AUTO-GENERATED SIMULATOR BY C2E2 */\n'

        # Include headers
        includes = ('# include <iostream>\n'
                    '# include <vector>\n'
                    '# include <boost/numeric/odeint.hpp>\n'
                    '# include <math.h>\n')

        # Set namespaces
        namespace = ('using namespace std;\n'
                    'using namespace boost::numeric::odeint;\n')

        # Create typedef
        typedef = 'typedef vector<double> state_t;\n'

        # Create integrator observer
        int_obs = ('//INTEGRATOR OBSERVER\n'
                'class IntObs {\n'
                '\tprivate:\n'
                '\t\tvector<state_t> &io_states;\n'
                '\t\tvector<double> &io_times;\n\n'
                '\tpublic:\n'
                '\t\tIntObs(vector<state_t> &states, vector<double> &times)\n'
                '\t\t\t: io_states(states), io_times(times) { }\n\n'
                '\t\tvoid operator()(const state_t &x, double t) {\n'
                '\t\t\tio_states.push_back(x);\n'
                '\t\t\tio_times.push_back(t);\n'
                '\t\t}\n'
                '};\n')

        # ODE functions
        odes = []
        for i, mode in enumerate(modes):
            ode = ['void ' + mode + '(const state_t &x, state_t &dxdt, const double t) {\n',
                dxdt[i] + '\n',
                '}\n']
            ode = ''.join(ode)
            odes.append(ode)
        odes = '//ODE FUNCTIONS\n' + '\n'.join(odes)

        # ODE function pointer
        ode_ptr = ('//ODE FUNCTION POINTER\n'
                'void (*rhs[' + str(len(modes)) + '])'
                '(const state_t &x, state_t &dxdt, const double t) =\n'
                '\t{' + ', '.join(modes) + '};\n')

        # Initialize variables
        init_vars = ['//VARIABLES',
                    'double ts, dt, te;',
                    'double abs_err, rel_err;',
                    'int cur_mode;',
                    'state_t x(' + str(len(vars)) + ');',
                    'vector<double> times;',
                    'vector<state_t> trace;']
        init_vars = '\t' + '\n\t'.join(init_vars) + '\n'

        # Read configuration file
        parse = ['//PARSING CONFIG',
                'cin >> ts;',
                'for (int i = 0; i < ' + str(len(vars)) + '; i++) {',
                '\tcin >> x[i];',
                '}',
                'cin >> abs_err >> rel_err >> dt >> te >> cur_mode;',
                'cur_mode--;']
        parse = '\t' + '\n\t'.join(parse) + '\n'

        # Integrate ODE
        integrate = ['//INTEGRATING']
        integrate.extend(integrator)
        integrate = '\t' + '\n\t'.join(integrate) + '\n'

        # FIXME have it only print the steps once without duplicate
        # Print step
        # FIXME code that I should be using
        """
        print_steps = ['//PRINTING TRACE',
                    'for (size_t i = 0; i <= steps; i++) {',
                    '\tcout << fixed;',
                    '\tcout << setprecision(9) << times[i];',
                    '\tfor (int j = 0; j < ' + str(len(vars)) + '; j++) {',
                    '\t\tcout << setprecision(10) << \' \' << trace[i][j];',
                    '\t}',
                    '\tcout << endl;',
                    '}']
        print_steps = '\t' + '\n\t'.join(print_steps)
        """

        # FIXME temporary to match CAPD simulator output
        print_steps = ['//PRINTING STEPS',
                    'for (size_t i = 0; i <= steps; i++) {',
                    '\tcout << fixed;',
                    '\tcout << setprecision(9) << times[i];',
                    '\tfor (int j = 0; j < ' + str(len(vars)) + '; j++) {',
                    '\t\tcout << setprecision(10) << \' \' << trace[i][j];',
                    '\t}',
                    '\tcout << endl;\n',
                    '\tif (i != 0 && i != steps) {',
                    '\t\tcout << fixed;',
                    '\t\tcout << setprecision(9) << times[i];',
                    '\t\tfor (int j = 0; j < ' + str(len(vars)) + '; j++) {',
                    '\t\t\tcout << setprecision(10) << \' \' << trace[i][j];',
                    '\t\t}',
                    '\t\tcout << endl;',
                    '\t}',
                    '}']
        print_steps = '\t' + '\n\t'.join(print_steps)
        
        # Generate main
        main = ['int main() {',
                init_vars,
                parse,
                integrate,
                print_steps,
                '}']
        main = '\n'.join(main)

        # Generate CPP file
        cpp_file = [auto_gen,
                    includes,
                    namespace,
                    typedef,
                    int_obs,
                    odes,
                    ode_ptr,
                    main]
        cpp_file = '\n'.join(cpp_file)

        # Write CPP file (NOTE: using absolute path here)
        file_path = Simulator.cur_path + file_path
        file = open(file_path, 'w')
        file.write(cpp_file)
        file.close()

    def compose_all(self):
        """ Compose Automaton """
        # TODO: Find a method to Parse Hybrid system(automata)
        # cls.parse(hybrid)

        # if not hybrid.parsed:
        #     Session.write("System not parsed. Exiting composition...\n")
        #     return

        # Session.write("Composing System...\n")
        
        automata_list = self.automata
        automata_list.reverse()
        while len(automata_list) > 1:
            automaton1 = automata_list.pop()
            automaton2 = automata_list.pop()
            automata_list.append(Simulator.compose(automaton1, automaton2))
        
        self.automata = automata_list
        # TODO: This!
        self.populateInvGuards()
        
        thinvarprop = ""
        thinvarlist = ""

        for automaton in self.automata:
            for var in automaton.local_var_names:
                if var in automaton.local_thinvar_names:
                    thinvarlist += var + "\n"
                    thinvarprop += "1\n"
                else:
                    thinvarprop += "0\n"

        writer = open(Simulator.cur_path + "../work-dir/ThinVarProp","w")
        writer.write(thinvarprop)
        writer.close()
        
        writer = open(Simulator.cur_path + "../work-dir/ThinVarList","w")
        writer.write(thinvarlist)
        writer.close()

        # NOTE: check if properties are parsed
        # cls.compose_properties(hybrid.automata[0], Session.hybrid.properties)
        
        # hybrid.parsed = False
        # cls.parse(hybrid)

        # Session.hybrid = hybrid
        # Session.hybrid.composed = True
        
        print("Composition complete.")
        
        return

    def populateInvGuards(self):
        """ Populate the guard and invariant dictionaries we use to generate PPL files """
        guardResets = defaultdict(list)
        invariants = defaultdict(list)
        # NOTE: NOT USED IN C2E2 
        # varList = self.local_var_names

        for automaton in self.automata:
            for m in automaton.mode_list:
                inv_eqs = []
                inv_vars = set()
                for inv in m.inv_list:
                    inv_eqs.append(inv.expr)
                    inv_vars = inv_vars.union(SymEq.vars_used(inv.expr))
                invariants[m.id] = (inv_eqs, inv_vars)

        for automaton in self.automata:
            for t in automaton.transition_list:
                g_eqs = t.guard.expr
                g_vars = SymEq.vars_used(t.guard.expr)
                action_eqs = []
                for act in t.actions:
                    action_eqs.extend(act.expr)
                guardResets[(t.source,t.destination)]\
                                                .append((g_eqs, action_eqs, g_vars))

        self.guardResets = dict(guardResets)
        self.invariants = dict(invariants)

    @staticmethod
    def compose(automaton1, automaton2):
        print("  Composing " + automaton1.name + " and " 
            + automaton2.name + "\n")
        composed = Automaton(automaton1.name + "_" + automaton2.name)
        # NOTE: NOT USED IN C2E2
        # m1_len = len(automaton1.modes)
        # m2_len = len(automaton2.modes)

        # modeid = 0;
        # Construct Cartesian product of modes
        for m1 in automaton1.modes:
            for m2 in automaton2.modes:
                m_name = m1.name + '_' + m2.name
                m_id = m1.id * automaton2.next_mode_id + m2.id
                # m_id=modeid
                # modeid = modeid+1
                m_initial = m1.initial and m2.initial
                cross_mode = Mode(name=m_name, isInitial=m_initial)

                # Replace input variables with corresponding output variables
                dai_dict = {}
                Simulator.construct_output_dict(cross_mode, m1.dais, dai_dict)
                Simulator.construct_output_dict(cross_mode, m2.dais, dai_dict)
                Simulator.replace_dais(cross_mode, dai_dict)

                # Resulting invariant is a conjunction
                for inv1 in m1.invariants:
                    cross_mode.add_invariant(inv1)
                for inv2 in m2.invariants:
                    cross_mode.add_invariant(inv2)

                # Add Mode with specified mode ID
                composed.add_mode(cross_mode, m_id)

        # Construct guards of composed automata. (a,b)->(a',b) & (a,b)->(a,b')
        # for all transitions a->a' and b->b'. Note there is no (a,b)->(a',b')
        trans_id = 0
        for t1 in automaton1.transitions:
            t_guard = t1.guard
            t_actions = t1.actions
          
            for m2 in automaton2.modes:
                i = m2.id
                t_src = t1.source * automaton2.next_mode_id + i
                t_dest = t1.destination * automaton2.next_mode_id + i
                cross_trans = Transition(guard=t_guard, 
                                         actions=t_actions,
                                         source=t_src, 
                                         destination=t_dest, 
                                         id=trans_id)
                composed.add_transition(cross_trans)
                trans_id += 1 

        for t2 in automaton2.transitions:
            t_guard = t2.guard
            t_actions = t2.actions
         
            for m1 in automaton1.modes:
                i = m1.id
                t_src = i * automaton2.next_mode_id + t2.source
                t_dest = i * automaton2.next_mode_id + t2.destination
                cross_trans = Transition(guard=t_guard, 
                                         actions=t_actions, 
                                         source=t_src, 
                                         destination=t_dest, 
                                         id=trans_id)
                composed.add_transition(cross_trans)
                trans_id += 1

        # Construct composed variables
        # composed.variables.local = automaton1.variables.local \
        #                            + automaton2.variables.local
        # composed.variables.output = automaton1.variables.output \
        #                             + automaton2.variables.output
        # composed.variables.input = automaton1.variables.input \
        #                            + automaton2.variables.input
        # composed.variables.input = [var for var in composed.variables.input \
        #                             if var not in composed.variables.output]
        composed.variable_list = automaton1.variable_list + automaton2.variable_list

        # Required to update the parents of the variable objects
        # composed.variables.update_parents(composed)

        print("  " + composed.name + " composition complete.\n")
        return composed

    @staticmethod
    def construct_output_dict(mode, dais, dai_dict):
        """
        Helper Function for compose:
        Construct dictionary mapping output variables to their RHS expressions
        """
        for dai in dais:
            lhs = str(dai.expr.lhs)
            rhs = str(dai.expr.rhs)
            if not lhs.endswith('_dot'):
                dai_dict[lhs] = rhs
            mode.add_dai(dai)

    @staticmethod
    def replace_dais(mode, dai_dict):
        """
        Helper Function for compose:
        Replace input variables with corresponding output variables
        """
        for dai in mode.dais:
            lhs = str(dai.expr.lhs)
            if lhs.endswith('_dot'):
                free_syms = dai.expr.rhs.free_symbols
                for var in free_syms:
                    var_name = str(var)
                    if var_name in dai_dict:
                        dai.expr = dai.expr.func(dai.expr.lhs, 
                            dai.expr.rhs.subs(var, dai_dict[var_name]))

                dai.raw = str(dai.expr.lhs) + ' = ' + str(dai.expr.rhs)

    # Generating Invariant & Guards Checking Funcitons
    def printHybridSimGuardsInvariants(self):
        self.printGuardsInvariants(Simulator.cur_path + "../work-dir/hybridSimGI.cpp", True)

    def printBloatedSimGuardsInvariants(self):
        self.printGuardsInvariants(Simulator.cur_path + "../work-dir/bloatedSimGI.cpp", False)

    def printGuardsInvariants(self, file_name, is_hybrid):
        checkFile = open(file_name, "w")
        codeString ="#include <ppl.hh>\n"
        codeString+="#include <iostream>\n"
        codeString+="#include <utility>\n"
        codeString+="#include <vector>\n"
        codeString+="#include <fstream>\n"
        codeString+="#include <typeinfo>\n\n"
        codeString+="using namespace std;\n\n"
        codeString+=self.printPoly()  # Constant
        codeString+=self.getMultFactorPt() # Local variables
        codeString+=self.getMultFactor() # Constant
        checkFile.write(codeString)
        checkFile.close()

        self.printInvariants(file_name, is_hybrid)
        self.printGuardResets(file_name, is_hybrid)

    def printInvariants(self, file_name, is_hybrid):
        """ 
        Generate PPL code to check invariants. 
        Only check against those variables in the invariant 
        """
        checkFile = open(file_name, "a")
        codeString="extern \"C\" bool invariantSatisfied(int curMode, double *ptLower, double *ptUpper){\n"
        codeString+="  NNC_Polyhedron box_poly;\n"
        codeString+="  double mult_factor = getMultFactor(" + ("ptUpper" if is_hybrid else "ptLower, ptUpper") + ");\n"
        for mode in self.invariants:
            eqs, varsUsed = self.invariants[mode]
            if not eqs:
                continue
            codeString+="  if(curMode=="+str(mode+1)+"){\n"
            codeString+= self.constructBoxHelper(varsUsed, "", 4, is_hybrid)
            codeString+="    Pointset_Powerset<NNC_Polyhedron> box(box_poly);\n"
            codeString+="    Pointset_Powerset<NNC_Polyhedron> invariant("+str(len(varsUsed))+",UNIVERSE);\n"
            codeString+="    Pointset_Powerset<NNC_Polyhedron> curInv;\n"
            codeString+="    NNC_Polyhedron curPoly;\n"
            codeString+="    Constraint_System cs;\n"
            for eq in eqs:
                codeString+="    curInv = Pointset_Powerset<NNC_Polyhedron>("+str(len(varsUsed))+",EMPTY);\n\n"
                for i,disjunct in enumerate(eq):
                    codeString+="    cs.set_space_dimension("+str(len(varsUsed))+");\n"
                    codeString+="    cs.insert("+str(disjunct)+");\n"
                    codeString+="    curPoly = NNC_Polyhedron(cs);\n"
                    codeString+="    curInv.add_disjunct(curPoly);\n"
                    codeString+="    cs.clear();\n\n"
                codeString+="    invariant.intersection_assign(curInv);\n\n"
            if is_hybrid: codeString+="    return invariant.contains(box);\n"
            else: codeString+="    return !(invariant.is_disjoint_from(box));\n"
            codeString+="  }\n"
        codeString+="  return true;\n"
        codeString+="}\n\n"
        checkFile.write(codeString)
        checkFile.close()

    def printGuardResets(self, file_name, is_hybrid):
        """ 
        Generates PPL code to check guards and create transitions.
        Only check against those variables in the guard. If reset exists, create 2*n dimensional space. 
        """
        varList = ["Simu_time"] + self.local_var_names
        resetVarList = [var+"_new" for var in varList]
        tempVarList = [var+"_temp" for var in varList]
        checkFile = open(file_name, "a")
        codeString="extern \"C\" vector<pair<NNC_Polyhedron, int> > hitsGuard(int curMode, double *ptLower, double *ptUpper){\n"
        codeString+="  vector<pair<NNC_Polyhedron, int> > toRet;\n"
        codeString+="  NNC_Polyhedron box_poly;\n"
        codeString+="  double mult_factor = getMultFactor(" + ("ptUpper" if is_hybrid else "ptLower, ptUpper") + ");\n"
        
        for key in self.guardResets:
            init = str(key[0]+1)
            dest = str(key[1]+1)
            for a,b,varsUsed in self.guardResets[key]:
                codeString+="  if(curMode=="+init+"){\n"
                codeString+= self.constructBoxHelper(varsUsed, "", 4, is_hybrid)
                thin_prop = -1
                for var in varsUsed:
                    if var in self.local_thinvar_names:
                        if thin_prop == 0:
                            thin_prop = -1
                            break
                        else:
                            thin_prop = 1
                    else:
                        if thin_prop == 1:
                            thin_prop = -1
                            break
                        else:
                            thin_prop = 0


                codeString+="    Constraint_System cs;\n"
                codeString+="    cs.set_space_dimension("+str(len(varsUsed))+");\n"    
    
                for guard_eq in a:
                    codeString+="    cs.insert("+str(guard_eq)+");\n"
                    
                codeString+="    NNC_Polyhedron guard(cs);\n"
                if is_hybrid: codeString+="    if(guard.contains(box_poly)){\n"
                else: codeString+="    if(!guard.is_disjoint_from(box_poly)){\n"
                codeString+= self.constructBoxHelper(tempVarList, "_temp", 6, is_hybrid)

                if is_hybrid and not b:
                    #codeString+="      toRet.push_back(make_pair(make_pair(box_poly,"+dest+"),"+str(thin_prop)+"));\n"
                    codeString+="      toRet.push_back(make_pair(box_poly,"+dest+"));\n"
                    codeString+="    }\n"
                    codeString+="  }\n"
                    continue

                codeString+="      Constraint_System cs_gd;\n"
                space_dim = 2*len(varList) if b else len(varList)
                codeString+="      cs_gd.set_space_dimension("+str(space_dim)+");\n"
                
                if b:
                    remVars = set(varList)
                    for i,var in enumerate(resetVarList):
                        codeString+="      Variable "+var+"("+str(len(varList)+i)+");\n"
                    codeString+="      box_poly.add_space_dimensions_and_embed("+str(len(varList))+");\n"
                    for reset_eq in b:
                        lhs, rhs = reset_eq.lhs, reset_eq.rhs
                        free_vars = list(lhs.free_symbols)
                        assert len(free_vars)==1
                        var = str(free_vars[0])
                        lhs = lhs.subs(var, var+'_new')
                        for v in self.local_var_names:
                            rhs = rhs.subs(v, v+'_temp')
                        codeString+="      cs_gd.insert("+str(lhs)+"=="+str(rhs)+");\n"
                        remVars.discard(var)
                    for var in remVars:
                        codeString+="      cs_gd.insert("+var+'_new'+'=='+var+"_temp);\n"
                
                if not is_hybrid:
                    for guard_eq in a:
                        for v in self.local_var_names:
                            guard_eq = guard_eq.subs(v, v+'_temp')
                        codeString+="      cs_gd.insert("+str(guard_eq)+");\n"
                codeString+="      NNC_Polyhedron guard_reset(cs_gd);\n"
                codeString+="      guard_reset.intersection_assign(box_poly);\n"

                if b:
                    codeString+="      Variables_Set vars;\n"
                    for var in varList:
                        codeString+="      vars.insert("+var+"_temp);\n"
                    codeString+="      guard_reset.remove_space_dimensions(vars);\n"
                
                #codeString+="      toRet.push_back(make_pair(make_pair(guard_reset,"+dest+"),"+str(thin_prop)+"));\n"
                codeString+="      toRet.push_back(make_pair(guard_reset,"+dest+"));\n"
                codeString+="    }\n"
                codeString+="  }\n"

        codeString+="  return toRet;\n"
        codeString+="}\n\n"
        checkFile.write(codeString)
        checkFile.close()

    # NOTE: NOT USED IN C2E2 ==============
    def constructBox(self, varList=None):
        if not varList:
            varList=["Simu_time"] + self.local_var_names
        codeString="NNC_Polyhedron constructBox(double *ptLower, double *ptUpper){\n"
        codeString+="  NNC_Polyhedron box_poly;\n"
        codeString+="  double mult_factor = getMultFactor(ptLower, ptUpper);\n"
        codeString+=self.constructBoxHelper(varList, "", 2)
        codeString+="  return box_poly;\n"
        codeString+="}\n\n"
        return codeString
    # =====================================

    def constructBoxHelper(self, varList, suffix, indent, is_point):
        indentation = " "*indent
        allVars = ["Simu_time"] + self.local_var_names

        codeString=indentation+"Constraint_System cs_box;\n"        
        for i,var in enumerate(varList):
            codeString+=indentation+"Variable "+var+"("+str(i)+");\n"
        codeString+="\n"

        for i,v in enumerate(allVars):
            if v+suffix in varList:
                var = v+suffix
                if is_point:
                    codeString+=indentation+"cs_box.insert(mult_factor*"+var+"==mult_factor*ptUpper["+str(i)+"]);\n"
                else:
                    codeString+=indentation+"if(ptLower["+str(i)+"]<ptUpper["+str(i)+"]){\n"
                    codeString+=indentation+"  cs_box.insert(mult_factor*"+var+">=mult_factor*ptLower["+str(i)+"]);\n"
                    codeString+=indentation+"  cs_box.insert(mult_factor*"+var+"<=mult_factor*ptUpper["+str(i)+"]);\n"
                    codeString+=indentation+"}\n"
                    codeString+=indentation+"else{\n"
                    codeString+=indentation+"  cs_box.insert(mult_factor*"+var+"<=mult_factor*ptLower["+str(i)+"]);\n"
                    codeString+=indentation+"  cs_box.insert(mult_factor*"+var+">=mult_factor*ptUpper["+str(i)+"]);\n"
                    codeString+=indentation+"}\n\n"
        codeString+=indentation+"box_poly = NNC_Polyhedron(cs_box);\n"
        return codeString

    def getMultFactorPt(self):
        """ Get factor to multiply doubles by because PPL only takes integers. """
        codeString="double getMultFactor(double *pt){\n"
        codeString+="  int multiplier=0, tmp_mul, str_len;\n"
        codeString+="  char buffer[100];\n"
        codeString+="  char *dot_loc;\n\n"
        mulLoop="  for(int i=0; i<"+str(len(self.local_var_names)+1)+"; i++){\n"
        mulLoop+="    sprintf(buffer, \"%lf\", pt[i]);\n"
        mulLoop+="    str_len = strlen(buffer);\n"
        mulLoop+="    dot_loc = strchr(buffer,'.');\n"
        mulLoop+="    if(dot_loc){\n"
        mulLoop+="      tmp_mul = (str_len-1)-(dot_loc-buffer);\n"
        mulLoop+="      if(tmp_mul>multiplier){\n"
        mulLoop+="        multiplier=tmp_mul;\n"
        mulLoop+="      }\n"
        mulLoop+="    }\n"
        mulLoop+="  }\n\n"
        codeString+=mulLoop
        codeString+="  return pow(10, multiplier);\n"
        codeString+="}\n\n"
        return codeString

    def getMultFactor(self):
        codeString="double getMultFactor(double *ptLower, double *ptUpper){\n"
        codeString+="  int lowerMult = getMultFactor(ptLower);\n"
        codeString+="  int upperMult = getMultFactor(ptUpper);\n"
        codeString+="  int multiplier = lowerMult > upperMult ? lowerMult : upperMult;\n"
        codeString+="  return multiplier;\n"
        codeString+="}\n\n"
        return codeString

    def printPoly(self):
        codeString="void print_box(NNC_Polyhedron poly){\n"
        codeString+="  Generator_System gs=poly.minimized_generators();\n"
        codeString+="  Generator_System::const_iterator i;\n"
        codeString+="  double divisor, dividend;\n"
        codeString+="  int dim;\n"
        codeString+="  cout << \"POLY: \" << endl;\n"
        codeString+="  for(i=gs.begin();i!=gs.end();++i){\n"
        codeString+="    if(i->is_point()){\n"
        codeString+="      divisor=mpz_get_d(i->divisor().get_mpz_t());\n"
        codeString+="      dim=int(i->space_dimension());\n"
        codeString+="      cout << \"POINT: \";\n"
        codeString+="      for(int j=0;j<dim;j++){\n"
        codeString+="        dividend=mpz_get_d(i->coefficient(Variable(j)).get_mpz_t());\n"
        codeString+="        cout<<dividend/divisor<<\" \";\n"
        codeString+="      }\n"
        codeString+="      cout<<endl;\n"
        codeString+="    }\n"
        codeString+="  }\n"
        codeString+="  cout << endl;\n"
        codeString+="}\n\n"
        return codeString


    # PROPERTIES
    @property
    def automata_names(self):
        """ Return all automaton's name in a list """
        if self.automata is not None:
            return [auto.name for auto in self.automata]
        else:
            return "No Automaton is assigned to this simulator."

    @property
    def local_var_names(self):
        var_names = []
        for automaton in self.automata:
            for name in automaton.local_var_names:
                var_names.append(name)
        return var_names

    @property
    def local_thinvar_names(self):
        thinvar_names = []
        for automaton in self.automata:
            for name in automaton.local_thinvar_names:
                thinvar_names.append(name)
        return thinvar_names

    @property
    def mode_names(self):
        """ Mode names of all automata """
        mode_names = []
        for automaton in self.automata:
            for mode_name in automaton.mode_names:
                mode_names.append(mode_name)
        return mode_names

    @property
    def modes(self):
        """ Mode() of all automata """
        modes = []
        for automaton in self.automata:
            for mode in automaton.modes:
                modes.append(mode)
        return modes     

    # CPP MODEL & COMPILE
    def initialize_cpp_model(self, config, sim_bool):

        model = Model()

        # Initialize set variables
        initial_mode = self.initial_set.initial_mode
        initial_eqns = self.initial_set.eqMatrix
        initial_matrix = SymEq.extract_matrix(self.initial_set.aMatrix, initial_eqns)
        initial_b = SymEq.extract_matrix(self.initial_set.bMatrix, initial_eqns)
        mode_names = self.mode_names
        initial_mode_idx = mode_names.index(initial_mode) + 1

        # Unsafe set variables
        unsafe_eqns = self.unsafe_set.eqMatrix
        unsafe_matrix = SymEq.extract_matrix(self.unsafe_set.aMatrix, unsafe_eqns)
        unsafe_b = SymEq.extract_matrix(self.unsafe_set.bMatrix, unsafe_eqns)

        # FIXME remove file readind and store in memory instead
        mode_linear = []
        gammas = []
        k_consts = []
        for m_i, m in enumerate(self.modes):
            fn = '../work-dir/jacobiannature' + str(m_i+1) + '.txt'
            fid = open(fn, 'r').read().split('\n')
            num_var = int(fid[0])

            if num_var == 0:
                m.linear = False

            if m.linear:
                list_var = []
                for i in range(num_var):
                    list_var.append(fid[i+1])

                eqn_pos = num_var + 1
                num_eqn = int(fid[eqn_pos])
                eqn_pos += 1

                list_eqn = []
                for i in range(num_eqn):
                    list_eqn.append(fid[eqn_pos+i])

                # FIXME see if we can avoid create functions dynamically

                codestring  = "def jcalc("
                codestring += "listofvar, "
                codestring += "listvalue"
                codestring += '):\n'
                codestring += " for i in range (len(listofvar)):\n"
                codestring += "   temp = listofvar[i]\n"
                codestring += "   rightside = '='+str(listvalue[i])\n"
                codestring += "   exec(temp+rightside)\n"
                codestring += " ret = []\n"
                for i in range (num_eqn):
                    codestring += " "
                    codestring += list_eqn[i]
                    codestring += '\n'
                    codestring += ' ret.append(Entry)\n'
                codestring += ' return ret'
                exec(codestring, globals())

                constant_jacobian = jcalc(list_var, np.ones((1, num_var))[0])
                constant_jacobian = np.reshape(constant_jacobian, (num_var, num_var))

                gamma_rate = np.linalg.eigvals(constant_jacobian).real
                gamma = max(gamma_rate)
                if(abs(max(gamma_rate)) < 0.00001):
                    gamma = 0
                k = np.linalg.norm(constant_jacobian)

            else:
                gamma = 0
                if self.property is not None:
                    k = self.property.k_value
                else:
                    # TODO: k=2000 for now, default value
                    k = 2000

            # Append calculated value
            mode_linear.append(int(m.linear))
            gammas.append(gamma)
            k_consts.append(k)

        # Unsigned integers
        model.set_dimensions(len(self.local_var_names))
        model.set_num_modes(len(self.modes))
        model.set_initial_mode_idx(initial_mode_idx)
        model.set_num_initial_eqns(len(initial_eqns))
        model.set_num_unsafe_eqns(len(unsafe_eqns))
        model.set_annot_type(3)

        # Integers
        if(self.refine_strat == DEF_STRAT):
            model.set_refine_strat(0)
        else:
            model.set_refine_strat(1)
            try:
                line = config["refine_order"]
                line = line.split(",")
                refine_order = [int(i) for i in line]
                # Int Vector
                refine_order_dv = IntegerVector()
                refine_order_dv[:] = refine_order
                model.set_refine_order(refine_order_dv)
            except:
                print("Refine order not provided, using default strategy")
                model.set_refine_strat(0)
                self.refine_strat = DEF_STRAT

        model.set_simulation_bool(sim_bool)

        # Strategy for blowing reachtube for Linear Model 
        if self.linear_strat == GSTAR:
            model.set_linear_strat(0)
        else:
            model.set_linear_strat(1)

        # Integer vectors
        mode_linear_dv  = IntegerVector()
        mode_linear_dv[:] = mode_linear
        model.set_mode_linear(mode_linear_dv)

        # Doubles
        model.set_abs_err(config['absolute_error'])
        model.set_rel_err(config['relative_error'])
        model.set_time_step(self.time_step)
        model.set_time_horizon(self.time_horizon)


        # Double vectors
        gammas_dv = DoubleVector()
        gammas_dv[:] = gammas
        k_consts_dv = DoubleVector()
        k_consts_dv[:] = k_consts
        initial_matrix_dv = DoubleVector()
        initial_matrix_dv[:] = initial_matrix
        initial_b_dv = DoubleVector()
        initial_b_dv[:] = initial_b
        unsafe_matrix_dv = DoubleVector()
        unsafe_matrix_dv[:] = unsafe_matrix
        unsafe_b_dv = DoubleVector()
        unsafe_b_dv[:] = unsafe_b

        model.set_gammas(gammas_dv)
        model.set_k_consts(k_consts_dv)
        model.set_initial_matrix(initial_matrix_dv)
        model.set_initial_b(initial_b_dv)
        model.set_unsafe_matrix(unsafe_matrix_dv)
        model.set_unsafe_b(unsafe_b_dv)

        model.set_visualize_filename('../work-dir/' + self.name)
        model.set_executable('../work-dir/' + config['simulator_name'])

        code_file=open("../work-dir/main.cpp","w+")
        code_file.write("#include <vector>\n")
        code_file.write('#include "model.hpp"\n')
        code_file.write("#define SIMU 1\n")
        code_file.write("#define VERI 0\n")
        code_file.write("using namespace std;")
        code_file.write("int main(){\n")
        code_file.write("    Model cpp_Model;\n")
        code_file.write("    cpp_Model.setDimensions("+str(len(self.local_var_names))+");\n")
        code_file.write("    cpp_Model.setNumModes("+str(len(self.modes))+");\n")
        code_file.write("    cpp_Model.setInitialModeIdx("+str(initial_mode_idx)+");\n")
        code_file.write("    cpp_Model.setNumInitialEqns("+str(len(initial_eqns))+");\n")
        code_file.write("    cpp_Model.setNumUnsafeEqns("+str(len(unsafe_eqns))+");\n")
        code_file.write("    cpp_Model.setAnnotType(3);\n")

        if(self.refine_strat == DEF_STRAT):
            code_file.write("    cpp_Model.setRefineStrat(0);\n")
        else:
            code_file.write("    cpp_Model.setRefineStrat(1);\n")
            code_file.write("    vector<int> refine_order_dv;\n")
            for i in range(len(refine_order)):
                code_file.write("    refine_order_dv.push_back("+ str(refine_order[i])+");\n")

        if self.linear_strat == GSTAR:
            code_file.write("    cpp_Model.setLinearStrat(0);\n")
        else:
            code_file.write("    cpp_Model.setLinearStrat(1);\n")

        code_file.write("    cpp_Model.setModeLinear(refine_order_dv);\n")

        code_file.write("    cpp_Model.setSimulationBool("+str(sim_bool)+");\n")

        code_file.write("    vector<int> mode_linear_dv;\n")
        for i in range(len(mode_linear)):
            code_file.write("    mode_linear_dv.push_back("+ str(mode_linear[i])+");\n")
        code_file.write("    cpp_Model.setModeLinear(mode_linear_dv);\n")

        code_file.write("    cpp_Model.setAbsError("+str(config['absolute_error'])+");\n")
        code_file.write("    cpp_Model.setRelError("+str(config['relative_error'])+");\n")
        code_file.write("    cpp_Model.setTimeStep("+str(self.time_step)+");\n")
        code_file.write("    cpp_Model.setTimeHorizon("+str(self.time_horizon)+");\n")

        code_file.write("    vector<double> gammas_dv;\n")
        for i in range(len(gammas_dv)):
            code_file.write("    gammas_dv.push_back("+ str(gammas[i])+");\n")
        code_file.write("    cpp_Model.setGammas(gammas_dv);\n")

        code_file.write("    vector<double> k_consts_dv;\n")
        for i in range(len(k_consts_dv)):
            code_file.write("    k_consts_dv.push_back("+ str(k_consts[i])+");\n")
        code_file.write("    cpp_Model.setKConsts(k_consts_dv);\n")

        code_file.write("    vector<double> initial_matrix_dv;\n")
        for i in range(len(initial_matrix_dv)):
            code_file.write("    initial_matrix_dv.push_back("+ str(initial_matrix[i])+");\n")
        code_file.write("    cpp_Model.setInitialMatrix(initial_matrix_dv);\n")

        code_file.write("    vector<double> initial_b_dv;\n")
        for i in range(len(initial_b_dv)):
            code_file.write("    initial_b_dv.push_back("+ str(initial_b[i])+");\n")
        code_file.write("    cpp_Model.setInitialB(initial_b_dv);\n")

        code_file.write("    vector<double> unsafe_matrix_dv;\n")
        for i in range(len(unsafe_matrix_dv)):
            code_file.write("    unsafe_matrix_dv.push_back("+ str(unsafe_matrix[i])+");\n")
        code_file.write("    cpp_Model.setUnsafeMatrix(unsafe_matrix_dv);\n")

        code_file.write("    vector<double> unsafe_b_dv;\n")
        for i in range(len(unsafe_b_dv)):
            code_file.write("    unsafe_b_dv.push_back("+ str(unsafe_b[i])+");\n")
        code_file.write("    cpp_Model.setUnsafeB(unsafe_b_dv);\n")

        code_file.write("    cpp_Model.setVisualizeFilename(\"../work-dir/" + self.name + "\");\n")
        code_file.write("    cpp_Model.setExecutable(\"../work-dir/" + config['simulator_name']+"\");\n")
        code_file.write("    cpp_Model.simulateVerify();\n")
        code_file.write("}\n\n")

        # Model() from CPP Backend
        return model

    def compile_executable(self, config):
        """ Compile Executable for Simulator """

        print("Compiling essential libraries for C2E2...")
        print("  Compilation may take a few minutes.")

        if((self.simulator == ODEINT_ADP) or (self.simulator == ODEINT_FIX)):
            print ("  Using ODEINT Simulator.\n")
            command_line = "g++ -w -O2 -std=c++11 simulator.cpp -o " + config['simulator_name']
            args = shlex.split(command_line)
            p = subprocess.Popen(args, cwd= "../work-dir")
            p.communicate()
        else:
            print("  Using CAPD Simulator.\n")
            command_line = "g++ -w -O2 simulator.cpp -o" + config['simulator_name'] + "`../capd/bin/capd-config --cflags --libs`"
            p = subprocess.Popen(command_line, cwd= "../work-dir", shell=True)
            p.communicate()

        command_line = "g++ -fPIC -shared hybridSimGI.cpp -o libhybridsim.so -lppl -lgmp"
        args = shlex.split(command_line)
        p = subprocess.Popen(args, cwd= "../work-dir")
        p.communicate()

        command_line = "g++ -fPIC -shared bloatedSimGI.cpp -o libbloatedsim.so -lppl -lgmp"
        args = shlex.split(command_line)
        p = subprocess.Popen(args, cwd= "../work-dir")
        p.communicate()

        command_line = "g++ -g -std=c++11 -fPIC -shared jaThin.cpp -o libJaThin.so -lppl -lgmp"
        args = shlex.split(command_line)
        p = subprocess.Popen(args, cwd= "../work-dir")
        p.communicate()
        
        print("Libraries successfully compiled.\n")
        return
