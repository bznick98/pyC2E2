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

from pyC2E2.constants import *
from pyC2E2.SymEq import SymEq

class Automaton:
    '''
    Data type represents a SINGLE automaton in hybrid model
    '''

    def __init__(self, name="default automaton"):
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
            "Variable List: {}\n".format(self.var_names()) + \
            "Mode List: {}\n".format(self.mode_names()) + \
            "Transition List: {}\n".format(self.transition_list)
        return DISPLAY_INFO

    # ADD to automaton
    def add_var(self, var):
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
        return

    def add_transition(self, tran):
        """ Add a transition to automaton """
        self.transition_list.append(tran)
        return

    # REMOVE from automaton


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
    
    # GETTERS: get automaton info
    def var_names(self):
        """ Return all variable names in a list """
        return [var.name for var in self.variable_list]
    
    def mode_names(self):
        """ Return all mode names in a list """
        return [mode.name for mode in self.mode_list]


class Variable:
    '''
    Automaton's Variable data type
    '''
    def __init__(self, name="default_variable",
                 type=REAL, scope='LOCAL_DATA'):
        self.name = name
        self.type = type
        self.scope = scope

    def __eq__(self, other):
        return (self.name==other.name\
                and self.type==other.type)


class ThinVariable:
    '''
    Automaton's Thin Variable data type
    '''
    def __init__(self, name="default_thinvariable", update_type="", 
                 type=REAL, scope='LOCAL_DATA'):
        self.name = name
        self.type = type
        self.scope = scope

    def __eq__(self, other):
        return (self.name==other.name\
                and self.type==other.type)

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

    # ADD to mode
    def add_dai(self, dai):
        '''
        Add a DAI equation to the mode
        '''
        self.dai_list.append(dai)
        return

    def add_invariant(self, inv):
        '''Add a Invariant to the mode'''
        self.inv_list.append(inv)
        return

    # REMOVE from mode


class Transition:
    '''
    guard   - node representing the guard for the transition
    actions - list of nodes representing the resets of the transition
    id      - unique identifier for the transition
    src     - the source of the transition
    dest    - the destination of the transition
    ''' 
    def __init__(self, guard, actions, id=-1, source=-1, destination=-1):
        # Call setters to ensure guard/action parents are set
        self.guard = guard
        self.actions = actions
        # Other variables, getters/setters not used
        self.id = id
        self.source = source
        self.destination = destination
        self.parent = None
    
    def __repr__(self):
        """ Display Transition in 'ModeId1 -> ModeId2' """
        return "{} -> {}".format(self.source, self.destination)


class DAI:
    '''
    Deterministic algebraic inequalities, a mode has DAIs
    '''

    def __init__(self, raw=None):
        
        self.raw = raw
        self.expr = None

class Invariant:
    '''
    Invariant data type, a mode has an invariant
    '''
    def __init__(self, raw=None):
        self.raw = raw
        self.expr = None
        self.parent = None

class Guard:
    '''
    Guard data type, a transition has a guard
    '''
    def __init__(self, raw=None):
        self.raw = raw
        self.expr = None
        self.parent = None

class Action:
    '''
    Action data type, a transition has an action
    '''
    def __init__(self, raw=None):

        self.raw = raw
        self.expr = None
        self.parent = None

    @property
    def name(self):
        return self.raw

    def parse(self):

        errors = []

        if self.raw is None:
            errors.append(('Action', self, "No Expression", None))
            return errors

        eqns = self.raw.split('&&')

        expr = []
        for eqn in eqns:
            constructed = SymEq.construct_eqn(eqn, True, True)
            if constructed is None:
                errors.append(('Action', self, "Invalid Expression", eqn))
            else:
                expr.append(constructed)
        self.expr = expr

        return errors