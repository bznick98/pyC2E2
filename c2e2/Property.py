# Property.py Class List:
# Property()

import re

from .constants import *
from .SymEq import SymEq

class Property():
    '''
    # TODO: Add Property() Definition
    '''
    def __init__(self):
        # Property values
        self.name = ''
        self.type = SAFETY
        self.time_step = 0.0
        self.time_horizon = 0.0
        self.k_value = 0.0
        self.initial_set_raw = ''
        self.initial_set = None
        self.unsafe_set_raw = ''
        self.unsafe_set = None
        self.simulator = ''

        # Property values valid flags (avoid repetitive checking?)
        self.name_valid = False
        self.type_valid = True
        self.time_step_valid = False
        self.time_horizon_valid = False
        self.k_value_valid = False
        self.initial_set_valid = False
        self.unsafe_set_valid = False

        # Misc
        self.status = "Not verified"
        self.result = ''
        self.is_visible = False

    def __repr__(self):
        """ Display Property Information """
        DISPLAY_INFO = \
            "\n=== Property() ===\n" + \
            "Name: {}           [Valid: {}]\n".format(self.name, self.name_valid) + \
            "Time Step: {}      [Valid: {}]\n".format(self.time_step, self.time_step_valid) + \
            "Time Horizon: {}   [Valid: {}]\n".format(self.time_horizon, self.time_horizon_valid) + \
            "K Value: {}        [Valid: {}]\n".format(self.k_value, self.k_value_valid) + \
            "Initial Set Raw: {} [Valid: {}]\n".format(self.initial_set_raw, self.initial_set_valid) + \
            "Unsafe Set Raw: {} [Valid: {}]\n".format(self.unsafe_set_raw, self.unsafe_set_valid) + \
            "Simulator: {}\n".format(self.simulator)
        return DISPLAY_INFO
            


    # Validate Property
    def is_valid(self):
        ''' Check all flags so you don't have to go through all validation again '''
        return (self.name_valid and\
                self.type_valid and\
                self.time_step_valid and\
                self.time_horizon_valid and\
                self.k_value_valid and\
                self.initial_set_valid and\
                self.unsafe_set_valid)

    def validate_initial(self):
        '''
        Only validate attributes that are available to validate during load property

        Args:
            None
        Returns:
            Set flags indicating the validity of each attribute
        '''

        self.validate_name()
        # self.validate_type(self.type)
        self.validate_time_step()
        self.validate_time_horizon()
        self.validate_k_value()
        return 

    def validate_all(self, automaton):
        '''
        validate ALL attributes in the Property class

        Args:
            automaton - the automaton that this property associated with,
                        because we need to validate is the initial mode exists 
                        in the automaton as well as variables in the initial set, 
                        we need to input this.
        
        Returns:
            Set flags indicating the validity of each attribute
        '''

        self.validate_name()
        # self.validate_type(self.type)
        self.validate_time_step()
        self.validate_time_horizon()
        self.validate_k_value()
        self.validate_initial_set(automaton)
        self.validate_unsafe_set(automaton)
        return 

    def validate_name(self):
        """ As long as the name is not Empty String """
        valid = True
        if self.name == '': valid = False
        # currently not checking repeated name  
        # set flag          
        self.name_valid = valid
        return valid
    
    def validate_time_step(self):
        """ valid = (is_number(time_step) and (time_step >= 0)) """
        valid = (self.time_step >= 0) and (is_number(self.time_step))
        # set flag
        self.time_step_valid = valid
        return valid

    def validate_time_horizon(self):
        """ valid = (is_number(horizon) and (horizon >= 0)) """
        valid = (self.time_horizon >= 0) and (is_number(self.time_horizon))
        # set flag
        self.time_horizon_valid = valid
        return valid
    
    def validate_k_value(self):
        """ valid = (is_number(k_value) and (k_value >= 0)) """
        valid = (self.k_value >= 0) and (is_number(self.k_value))
        # set flag
        self.k_value_valid = valid
        return valid

    def validate_initial_set(self, automaton):
        """
        Validating Initial Set (In String Format) of the Property 
        
        Args:
            automaton - the automaton that this property associated with,
                        because we need to validate is the initial mode exists 
                        in the automaton as well as variables in the initial set, 
                        we need to input this.
        Return
            If valid, set initial_set_valid flag to True
        """
        # Regex building blocks
        flt = '(-?(\d*\.?\d+))'
        int = '(-?(\d+))'
        term = '(' + flt + '|(' + int + '\s*/\s*' + int + '))'
        var = '([a-zA-Z]\w*)'
        mode = '(' + var + ':)'
        eql = '((<=?)|(>=?)|(==))'
        expr = '(' + term + '?\s*' + var + '(\s*[-\+]?\s*((' + term + '?\s*' + var + ')|' + term + '))*)'
        eqn = '(' + expr + '\s*' + eql + '\s*' + term + ')'
        eqns = '(' + eqn + '(\s*&&\s*' + eqn + ')*)'
        init_set = '(\s*' + mode + '\s*' + eqns + '\s*){1}$'
        unsafe_set = '(\s*' + eqns + '\s*)+$'

        # Regex strings
        re_mode = mode
        re_var = var
        re_is = init_set
        re_us = unsafe_set

        initial_set_raw = self.initial_set_raw

        match = re.match(re_is, initial_set_raw)

        if(match == None):
            self.initial_set = None
            self.initial_set_valid = False
            raise Exception("[Initial Set Error]: Incorrect Syntax")
        else:
            is_sep = initial_set_raw.split(':')

            # Validate Mode
            mode = re.search(re_var, is_sep[0]).group(0)
            mode_list = [mode.name for mode in automaton.mode_list]

            if(mode not in mode_list):
                self.initial_set = None
                self.initial_set_valid = False
                raise Exception("[Initial Set Error]: No matching modes")

            # Validate Vars
            vars_ = re.findall(re_var, is_sep[1])
            var_name_list = automaton.var_names()
            var_union= set(vars_) | set(var_name_list)

            if(len(var_union) > len(var_name_list)):
                self.initial_set = None
                self.initial_set_valid = False
                raise Exception("[Initial Set Error]: Variable mismatch")

            # Parse equations
            a_m, b_m, eq_m = SymEq.get_eqn_matrix(is_sep[1], var_name_list)
            bounded = SymEq.check_boundedness(a_m, b_m, eq_m, var_name_list)

            if(is_sep[1].count('>')!= is_sep[1].count('<')):
                bounded = False
            
            if bounded:
                self.initial_set = [is_sep[0], a_m, b_m, eq_m]
                self.initial_set_valid = True
                return
            else:
                self.initial_set = None
                self.initial_set_valid = False
                raise Exception("[Initial Set Error]: Set unbounded")


    def validate_unsafe_set(self, automaton):
        """
        Validating Unsafe Set
        Args:
            automaton - the automaton that this property associated with,
                        because we need to validate if variables match those,
                        inthe unsafe_set, so need to input this.
        Return
            If valid, set unsafe_set_valid flag to True
        """
        # Regex building blocks
        flt = '(-?(\d*\.?\d+))'
        int = '(-?(\d+))'
        term = '(' + flt + '|(' + int + '\s*/\s*' + int + '))'
        var = '([a-zA-Z]\w*)'
        mode = '(' + var + ':)'
        eql = '((<=?)|(>=?)|(==))'
        expr = '(' + term + '?\s*' + var + '(\s*[-\+]?\s*((' + term + '?\s*' + var + ')|' + term + '))*)'
        eqn = '(' + expr + '\s*' + eql + '\s*' + term + ')'
        eqns = '(' + eqn + '(\s*&&\s*' + eqn + ')*)'
        init_set = '(\s*' + mode + '\s*' + eqns + '\s*){1}$'
        unsafe_set_reg = '(\s*' + eqns + '\s*)+$'

        # Regex strings
        re_mode = mode
        re_var = var
        re_is = init_set
        re_us = unsafe_set_reg

        # Check if input is valid
        if SymEq.check_expression_syntax(self.unsafe_set_raw, "Unsafe Set"):
            self.unsafe_set = None
            self.unsafe_set_valid = False
            raise Exception("[Unsafe Set Error]: Incorrect Syntax")
        
        # Validate vars
        else:
            vars_ = re.findall(re_var, self.unsafe_set_raw)
            var_list = automaton.var_names()
            var_union = set(vars_) | set(var_list)

            if(len(var_union) > len(var_list)):
                self.unsafe_set = None
                self.unsafe_set_valid = False
                raise Exception("[Unsafe Set Error]: Variable mismatch")

            self.unsafe_set = SymEq.get_eqn_matrix(self.unsafe_set_raw, var_list)
            self.unsafe_set_valid = True
            return


# Helper Functions:
def is_number(number):
    """ Helper Function: Checking the argument is a Number(int float complex) """
    return isinstance(number, (int, float, complex))