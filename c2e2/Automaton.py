# Automaton.py Class List:
# Autonmation()
# Variable()
# ThinVariable()
# Mode()
# Transition()
# DAI()
# Invariant()
# Guard()
# Action()

import sympy
import warnings

from .constants import *
from .SymEq import SymEq

class Automaton:
    '''
    A single Automaton representation.

    An Automaton contains all necessary information such as the dynamics of the hybrid model,
    To do Simulation/Reachability Analysis, Information such as Initial Set and Time Step,
    Time Horizon, etc. need to be provided.

    name             - the name of the automaton
    variable_list    - all variables in the automaton, a list of Variable()
    thinvar_list     - all Thin-variables in the automaton, a list of ThinVariable()
    mode_list        - all modes of the automaton, a list of Mode()
    transistion_list - all transitions of the automaton, a Transition() is from old_mode to new_mode
    '''

    def __init__(self, name:str="default automaton"):
        self.name = name
        self.variable_list = []
        self.thinvar_list = []      # not too sure what thin var is
        self.mode_list = []
        self.transition_list = []

    def __repr__(self):
        """ Display Automaton Info """
        DISPLAY_INFO = \
            "\n=== Automaton() ===\n" + \
            "Name: {}\n".format(self.name) + \
            "Variable List: {}\n".format(self.var_names) + \
            "Mode List: {}\n".format(self.mode_names) + \
            "Transition List: {}\n".format(self.transition_list)
        return DISPLAY_INFO

    # Setters
    def add_variable(self, var):
        """ Add a variable to the automaton """
        self.variable_list.append(var)
        return
    
    def add_thinvar(self, thinvar):
        """ Add a thin variable to the automaton """
        self.thinvar_list.append(thinvar)
        return

    def add_mode(self, mode, id=None):
        """ Add a mode to automaton """
        self.mode_list.append(mode)
        # TODO: PARSE MODE & ALSO ADD VARIABLES TO Automaton
        return

    def add_transition(self, tran):
        """ Add a transition to automaton """
        self.transition_list.append(tran)
        return

        
    # Getters
    # Variables
    @property
    def vars(self):
        """ Return all variable in a list """
        return [var for var in self.variable_list]

    @property
    def var_names(self):
        """ Return all variable names in a list """
        return [var.name for var in self.variable_list]

    @property
    def local_vars(self):
        """ Return all variable whose scope is 'LOCAL' in a list """
        local_vars = []
        for var in self.variable_list:
            if var.scope == LOCAL:
                local_vars.append(var)
        return local_vars
    
    @property
    def local_var_names(self):
        """ Return all variable name whose scope is 'LOCAL' in a list """
        return [var.name for var in self.local_vars]

    @property
    def input_vars(self):
        """ Return all variable whose scope is 'INPUT' in a list """
        input_vars = []
        for var in self.variable_list:
            if var.scope == INPUT:
                input_vars.append(var)
        return input_vars
    
    @property
    def input_var_names(self):
        """ Return all variable name whose scope is 'INPUT' in a list """
        return [var.name for var in self.input_vars]

    @property
    def output_vars(self):
        """ Return all variable whose scope is 'OUTPUT' in a list """
        output_vars = []
        for var in self.variable_list:
            if var.scope == OUTPUT:
                output_vars.append(var)
        return output_vars
    
    @property
    def output_var_names(self):
        """ Return all variable name whose scope is 'OUTPUT' in a list """
        return [var.name for var in self.output_vars]

    # ThinVariables
    @property
    def local_thinvars(self):
        """ Return all ThinVariable whose scope is 'LOCAL' in a list """
        local_vars = []
        for var in self.thinvar_list:
            if var.scope == LOCAL:
                local_vars.append(var)
        return local_vars
    
    @property
    def local_thinvar_names(self):
        """ Return all variable name whose scope is 'LOCAL' in a list """
        return [var.name for var in self.local_thinvars]

    @property
    def input_thinvars(self):
        """ Return all ThinVariable whose scope is 'INPUT' in a list """
        input_vars = []
        for var in self.thinvar_list:
            if var.scope == INPUT:
                input_vars.append(var)
        return input_vars
    
    @property
    def input_thinvar_names(self):
        """ Return all variable name whose scope is 'INPUT' in a list """
        return [var.name for var in self.input_thinvars]

    @property
    def output_thinvars(self):
        """ Return all variable whose scope is 'OUTPUT' in a list """
        output_vars = []
        for var in self.thinvar_list:
            if var.scope == OUTPUT:
                output_vars.append(var)
        return output_vars
    
    @property
    def output_thinvar_names(self):
        """ Return all variable name whose scope is 'OUTPUT' in a list """
        return [var.name for var in self.output_thinvars]

    # Modes
    @property
    def modes(self):
        """ Return all modes in a list """
        return [mode for mode in self.mode_list]

    @property
    def mode_names(self):
        """ Return all mode names in a list """
        return [mode.name for mode in self.mode_list]

    # TODO: Functions to REMOVE attributes from automaton

    # CHECKING automton validity
    def verify_mode_ids(self):
        """ Verify mode ids are unique within an Automaton """
        id_set = set()
        for mode in self.mode_list:
            if mode.id in id_set:
                return False
            id_set.add(mode.id)
        return True

    def verify_mode_names(self):
        """ Verify mode names are unqiue with an Automaton """
        name_set = set()
        for mode in self.mode_list:
            name = mode.name.strip()
            if name in name_set:
                return False
            name_set.add(name)
        return True

    def verify_transition_src_dest(self):
        """ Verify transition src and destination IDs """
        id_set = set()
        for mode in self.mode_list:
            id_set.add(mode.id)

        for tran in self.transition_list:
            if ((tran.source not in id_set) or 
                (tran.destination not in id_set)):
                return False
        return True


class Variable:
    '''
    Represent a single Variable belongs to the Automaton.

    name    - name of the Variable
    type    - REAL/INTEGER, data type of the Variable
    scope   - LOCAL/INPUT/OUTPUT, scope of the Variable
    '''
    def __init__(self, name="default_variable",
                 type=REAL, scope=LOCAL):
        self.name = name
        self.type = type
        self.scope = scope

    def __eq__(self, other):
        return (self.name==other.name\
                and self.type==other.type)
    
    def __repr__(self):
        DISPLAY_INFO = "\n=== Variable() ===\n"
        DISPLAY_INFO += "Name: {}".format(self.name)
        DISPLAY_INFO += "Type: {}".format(self.type)
        DISPLAY_INFO += "Scope: {}".format(self.scope)


class ThinVariable:
    '''
    Represent a single ThinVariable belongs to the Automaton.

    name    - name of the ThinVariable
    type    - REAL/INTEGER, data type of the Variable
    scope   - LOCAL/INPUT/OUTPUT, scope of the Variable
    '''
    def __init__(self, name="default_thinvariable", update_type="", 
                 type=REAL, scope='LOCAL_DATA'):
        self.name = name
        self.type = type
        self.scope = scope

    def __eq__(self, other):
        return (self.name==other.name\
                and self.type==other.type)
    
    def __repr__(self):
        DISPLAY_INFO = "\n=== ThinVariable() ===\n"
        DISPLAY_INFO += "Name: {}".format(self.name)
        DISPLAY_INFO += "Type: {}".format(self.type)
        DISPLAY_INFO += "Scope: {}".format(self.scope)


class Mode:
    '''
    name        - name of the mode (can be different from the one specified in 
                  simulink/stateflow
    id          - uniquely identifies mode. if mode created through HyIR, then id 
                  corresponds to its index in the HyIR's modes list
    isInitial   - True if this is the starting mode of the hybrid automata, 
                  false otherwise
    inv_list    - list of Invariant objects representing the mode's invariants
    dai_list    - list of DAI objects representing the mode's governing differential 
                  equations
    '''
    def __init__(self, name="default_mode", id=-1, isInitial=False):
        # Call setters to ensure transition/dai parents are set
        self.inv_list = []
        self.dai_list = []
        # Other variables, getters/setters not used
        self.name = name
        self.id = id
        self.isInitial = isInitial
        # self.initialConditions = []
        self.linear = True
        self.parent = None

    def __repr__(self):
        DISPLAY_INFO = "\n=== Mode() ===\n"
        DISPLAY_INFO += "Mode Name: {}".format(self.name)
        if self.isInitial:
            DISPLAY_INFO += " (Initial Mode)\n"
        else:
            DISPLAY_INFO += "\n"
        DISPLAY_INFO += "Mode ID: {}\n".format(self.id)
        if self.linear is None:
            DISPLAY_INFO += "Linearity: Unknown\n"
        elif self.linear:
            DISPLAY_INFO += "Linearity: Linear\n"
        else:
            DISPLAY_INFO += "Linearity: Non-Linear\n"
        # DAIs
        DISPLAY_INFO += "DAI Equations: "
        if len(self.dai_list) == 0:
            DISPLAY_INFO += "(None)\n"
        else:
            DISPLAY_INFO += "\n"
            for dai in self.dai_list:
                DISPLAY_INFO += "  {}\n".format(dai.raw)
        # Invariants
        DISPLAY_INFO += "Invariants: "
        if len(self.inv_list) == 0:
            DISPLAY_INFO += "(None)\n"
        else:
            DISPLAY_INFO += "\n"
            for inv in self.inv_list:
                DISPLAY_INFO += "  {}\n".format(inv.raw)
        return DISPLAY_INFO

    @property
    def dais(self):
        """ Same as self.dai_list """
        return self.dai_list
        
    @property
    def invariants(self):
        """ Same as self.inv_list """
        return self.inv_list

    # @property
    # def vars(self):
    #     """ Return all vars in a list of Variable() """
    #     var_list = []
    #     # get all raw expressions
    #     raw_dai_list = []
    #     for dai in self.dai_list:
    #         raw_dai_list.append(dai.raw)
    #     # extract all variables(repeats)
    #     for raw_expr in raw_dai_list:
    #         expr_var_list = SymEq.get_var_list(raw_expr)
    #         for var in expr_var_list:
    #             var_list.append(var)
    #     # remove repeat
    #     var_list = list(set(var_list))
    #     # convert to Variable()
    #     for var in var_list:
    #         var = Variable(var)
    #     print(var_list)
    #     return var_list
        
    # ADD to mode
    def add_dai(self, dai):
        '''
        Add a DAI equation to the mode
        '''
        if not isinstance(dai, DAI):
            raise Exception("[Mode's DAI Error]: add_dai() only accept DAI object.")
        elif dai.raw is None:
            # warnings will not affect program running
            warnings.warn("[Mode's DAI Warning]: An empty DAI equation, ignored")
        else:
            self.dai_list.append(dai)
        return

    def add_invariant(self, inv):
        '''Add a Invariant to the mode'''
        if isinstance(inv, Invariant):
            self.inv_list.append(inv)
        else:
            raise Exception("[Mode's Invariant Error]: add_invariant() only accept Invariant object.")
        return

    # REMOVE from mode


    # PARSING and CHECKING
    def parse(self):
        """ Parse DAI equation and Invariant Equations """
        # DAI
        self.parse_all_dais()
        # Invariants
        self.parse_all_invariants()
        return

    def parse_all_dais(self):
        """ Parse DAIs """
        self.linear = True        
        for dai in self.dai_list:
            dai.parse()
            if (self.linear and (dai.expr is not None)):
                self.linear = SymEq.is_linear(dai.expr.rhs)

        # TODO: Parse DAI variables
        # self.parse_dai_vars()

    def parse_all_invariants(self):
        """ Parse Invariants """
        for inv in self.inv_list:
            p = inv.parse()
            errors += inv.parse()


    # def parse_dai_vars(self):
    #     """ 
    #     Parse variables used in DAI equations
    #     - Output variables on left-hand side must not end with _dot
    #     - Local variables on left-hand side must end with _dot
    #     - Input variables must not be used on the left-hand side
        
    #     Assumptions
    #     - <DAI>.expr contains only one symbol on the left-hand side
    #         - This check is done during <DAI>.parse()
    #     """

    #     errors = []
    #     local_vars = set(self.parent.local_var_names)
    #     output_vars = set(self.parent.output_var_names)
    #     input_vars = set(self.parent.input_var_names)

    #     for dai in self.dai_list:

    #         if dai.expr is None:
    #             continue

    #         lhs = str(dai.expr.lhs)
    #         var_ = lhs.replace('_dot', '')
            
    #         if var_ in local_vars:
    #             if not lhs.endswith('_dot'):
    #                 errors.append(('Flow', dai, "Incorrect Variable Usage",
    #                                "Local variable equations must end with "
    #                                + "_dot"))
    #             local_vars.remove(var_)
    #         elif var_ in output_vars:
    #             if lhs.endswith('_dot'):
    #                 errors.append(('Flow', dai, "Incorrect Variable Usage", 
    #                                "Output variable equations must not end "
    #                                + "with _dot"))
    #             output_vars.remove(var_)
    #         elif var_ in input_vars:
    #             errors.append(('Flow', dai, "Incorrect Variable Usage",
    #                            "Input variables should not appear on the "
    #                            "left-hand side of an equation"))

    #     for var_ in local_vars.union(output_vars):
    #         errors.append(('WARNING', self, "Unused variable: " + var_, None))
        
    #     return errors


class Transition:
    '''
    guard   - node representing the guard for the transition
    actions - list of nodes representing the resets of the transition
    id      - unique identifier for the transition
    src     - the source of the transition
    dest    - the destination of the transition
    parent  - TODO: Add definition for this 
    ''' 
    def __init__(self, guard, actions=[], id=0, src_mode=0, dst_mode=0):
        # Call setters to ensure guard/action parents are set
        self.guard = guard
        self.actions = actions
        # Other variables, getters/setters not used
        self.id = id
        self.source = src_mode
        self.destination = dst_mode
        self.parent = None
    
    def __repr__(self):
        """ Display Transition in 'ModeId1 -> ModeId2' """
        return "{} -> {}".format(self.source, self.destination)


class DAI:
    '''
    Deterministic algebraic inequalities, to represent differential equations.

    raw     - Raw (in string) Expressions
    expr    - Equations in the form of SymPy Type
    '''
    def __init__(self, raw=None):
        # Same as set_expression(raw)
        self.raw = raw
        self.expr = self.parse()

    def __call__(self, raw=None):
        """ Calling the object directly will set the input expression to DAI """
        self.set_expression(raw)

    def __repr__(self):
        DISPLAY_INFO = "\n=== DAI() ===\n"
        DISPLAY_INFO += "Raw Expression: {}\n".format(self.raw)
        DISPLAY_INFO += "SymPy Equation: {}\n".format(self.expr)
        return DISPLAY_INFO

    def set_expression(self, raw):
        """ call DAI object with a expression will set the expression to the DAI """
        self.raw = raw
        self.expr = self.parse()
    
    def parse(self):
        """ Parsing a raw expression and set its SymPy expression to the object """
        if self.raw is not None:
            # Construct SymPy Equation
            expr = SymEq.construct_eqn(self.raw, is_eq=True, rationalize=False)
            # Check lhs single symbol
            if type(expr.lhs) is not sympy.Symbol:
                self.raw = None
                raise Exception('Flow:' + self.raw + " -> Invalid Expession: Left hand side" +
                            " of equation must be a single symbol.")
            return expr
        else:
            return None


class Invariant:
    '''
    Invariant is a set of inequalities that has to be strictly satisfied during automaton's time period.

    raw     - Raw (in string) Expressions
    expr    - Equations in the form of SymPy Type
    parent  - TODO: Add definition for this
    '''
    def __init__(self, raw=None):
        self.raw = raw
        self.expr = None
        self.parent = None

        # parse raw expression
        self.parse()

    def __repr__(self):
        DISPLAY_INFO = "\n=== Invariant() ===\n"
        DISPLAY_INFO += "Raw Expression: {}\n".format(self.raw)
        DISPLAY_INFO += "SymPy Equation: {}\n".format(self.expr)
        return DISPLAY_INFO

    def __call__(self, raw):
        """ Update/Assign new expressions by directly calling an Invariant object. """
        self.raw = raw
        self.parse()

    def parse(self):
        """ Parsing Raw string and convert it to Expression """
        if self.raw is None:
            warnings.warn("[Invariant Warning]: No Expression.")
            
        eqns = self.raw.split('||')
        expr = []
        for eqn in eqns:
            constructed = SymEq.construct_eqn(eqn, False, True)
            if constructed is None:
                temp_raw = self.raw
                self.raw = None
                raise Exception("[Invariant Error]: Invalid Syntax for Invariant: {}".format(temp_raw))
            else:
                expr.append(constructed)
        self.expr = expr

        # Filter out equations that evaluate to False
        self.expr = filter(lambda eqn: eqn is not False, self.expr)
        self.expr = list(self.expr)
        if True in self.expr: 
            temp_raw = self.raw
            self.raw = None
            self.expr = None
            raise Exception("[Invariant Error]: Redundant Invariant: {}".format(temp_raw))

        return


class Guard:
    '''
    Guard is a set of inequalities that when satisfied (depending on non-deterministic / urgent transitions),
    the automaton is going to switch mode.

    raw     - Raw (in string) Expressions
    expr    - Equations in the form of SymPy Type
    parent  - TODO: Add definition for this
    '''
    def __init__(self, raw=None):
        self.raw = raw
        self.expr = None
        self.parent = None

        # parse raw expression
        self.parse()

    def __repr__(self):
        DISPLAY_INFO = "\n=== Guard() ===\n"
        DISPLAY_INFO += "Raw Expression: {}\n".format(self.raw)
        DISPLAY_INFO += "SymPy Equation: {}\n".format(self.expr)
        return DISPLAY_INFO

    def __call__(self, raw):
        """ Update/Assign new expressions by directly calling a Guard object. """
        self.raw = raw
        self.parse()

    def parse(self):
        if self.raw is None:
            warnings.warn("[Guard Warning]: No Expression.")

        eqns = self.raw.split('&&')
        expr = []
        for eqn in eqns:
            constructed = SymEq.construct_eqn(eqn, False, True)
            if constructed is None:
                temp_raw = self.raw
                self.raw = None
                raise Exception("[Guard Error]: Invalid Syntax for Guard: {}". format(temp_raw))
            else:
                expr.append(constructed)
        self.expr = expr

        # Filter out equations that evaluate to True
        self.expr = filter(lambda eqn: eqn is not True, self.expr)
        self.expr = list(self.expr)
        if False in self.expr: 
            temp_raw = self.raw
            self.raw = None
            self.expr = None
            raise Exception("[Guard Error]: Redundant Guard: {}".format(temp_raw))
        
        return


class Action:
    '''
    Action will assign values to automaton's variables if the guards were hit.

    raw     - Raw (in string) Expressions
    expr    - Equations in the form of SymPy Type
    parent  - TODO: Add definition for this
    '''
    def __init__(self, raw=None):

        self.raw = raw
        self.expr = None
        self.parent = None
        # parsing raw expression
        self.parse()
    
    def __repr__(self):
        DISPLAY_INFO = "\n=== Action() ===\n"
        DISPLAY_INFO += "Raw Expression: {}\n".format(self.raw)
        DISPLAY_INFO += "SymPy Equation: {}\n".format(self.expr)
        return DISPLAY_INFO

    @property
    def name(self):
        return self.raw

    def parse(self):
        if self.raw is None:
            raise Exception("[Action Error]: No Expression found.")

        eqns = self.raw.split('&&')

        expr = []
        for eqn in eqns:
            constructed = SymEq.construct_eqn(eqn, True, True)
            if constructed is None:
                raise Exception("[Action Error]: Invalid Expression")
            else:
                expr.append(constructed)
        self.expr = expr
        return