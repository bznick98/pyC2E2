# Set.py Class List:
# Set() - An Abstract Class
# RectangleSet()

import re

from abc import ABC, abstractmethod
from .SymEq import SymEq

class Set(ABC):
    '''
    Set representation, an abstract class
    '''
    def __init__(self):
        self.name = "This is an abstract Set"

    @abstractmethod
    def raw_to_matrix(self):
        pass

class RectangleSet(Set):
    '''
    Rectangle Set:
    Example: Suppose there is an Initial Set raw expression:
    x >= 14 && x<=14.95 && y >= 3.25 && y<=3.3 && z == 0

    aMatrix will be: 
    [   [10, 0  , 0], [100, 0 , 0],     - for x
        [0 , 100, 0], [0  , 10, 0],     - for y
        [0 , 0  , 1], [0  , 0 , 1]  ]   - for z
             |              |
        for >/>=        for </<=

    bMatrix will be: 
    [   [140],      [1495],      - for x
        [325],      [33],        - for y
        [0],        [0]     ]    - for z
         |           |
       for >/>=   for </<=

    eqMatrix will be:
    [   ['>='],     ['<='],     - for x
        ['>='],     ['<='],     - for y
        ['>='],     ['<=']  ]   - for z
          |           |
        for >/>=   for </<=
    '''
    def __init__(self, expression=None, Eq_Matrix=None):
        """ Construct Set Using either expressions or Eq Matrix(a, b, eq) """
        self.name = "RectangleSet Default Name"
        self.raw_expression = expression
        self.initial_mode = None
        self.aMatrix = None
        self.bMatrix = None
        self.eqMatrix = None

        # process input expressions if user gives one
        if expression is not None:
            self.set_expression(expression)
        elif Eq_Matrix is not None:
            self.aMatrix = Eq_Matrix[0]
            self.bMatrix = Eq_Matrix[1]
            self.eqMatrix = Eq_Matrix[2]
    
    def __call__(self, expression):
        """ Directly Calling the object will set input expression """
        self.set_expression(expression)

    def __repr__(self):
        DISPLAY_INFO = '\n=== RectangleSet() ===\n'
        if self.raw_expression is not None:
            DISPLAY_INFO += '[Raw Expression]: '
            DISPLAY_INFO += self.raw_expression + "\n"
            if self.initial_mode is not None:
                DISPLAY_INFO += "[Initial Mode]: "
                DISPLAY_INFO += self.initial_mode +"\n"
        else:
            DISPLAY_INFO += "[Empty Set]: No Raw Expressions Available."
        return DISPLAY_INFO

    def set_expression(self, expression):
        """ Add expressions to the Rectangle Set and convert them to matrix representations """

        self.raw_expression = expression

        # Check if expression contains initial mode
        splited = self.raw_expression.split(":")
        if len(splited) == 1:
            # CHECKING new expression validity
            self.check_syntax(isInitialSet=False)
            # string contains initial mode
            self.initial_mode = None
            # no initial mode is provided(or maybe it's unsafe set, doesn't need a initial mode)
            self.aMatrix, self.bMatrix, self.eqMatrix = self.raw_to_matrix()
        else:
            # CHECKING new expression validity
            self.check_syntax(isInitialSet=True)
            # string contains initial mode
            self.initial_mode = splited[0].strip()
            self.raw_expression = splited[1].strip()
            # convert the expression to Matrix
            self.aMatrix, self.bMatrix, self.eqMatrix = self.raw_to_matrix()
        
        if self.initial_mode is not None:
        # if expression has mode information, then it must be initial set, so boundedness has to be checked
            self.check_boundedness()

    def raw_to_matrix(self):
        """ Using SymEq's get_eqn_matrix function """
        # extract variable list (initial set guarantees to include all vars)
        expression = self.raw_expression
        varList = SymEq.get_var_list(expression)
        return SymEq.get_eqn_matrix(expression, varList)


    # ERROR CHECKING: Syntax
    def check_syntax(self, isInitialSet=False):
        """ Using SymEq's Check Expression Syntax """
        expression = self.raw_expression
        if SymEq.check_expression_syntax(expression, isInitialSet):
            return True
        else:
            self.raw_expressions = None
            raise Exception("[RectangleSet ERROR]: Syntax NOT valid.")

    # ERROR CHECKING: Boundedness
    def check_boundedness(self):
        """ Using SymEq's Check boundedness, Mainly for Initial Set """
        if SymEq.check_boundedness(self.aMatrix,
                                    self.bMatrix,
                                    self.eqMatrix,
                                    SymEq.get_var_list(self.raw_expression)):
            return True
        else:
            raise Exception("[RectangleSet ERROR]: (Initial) Set NOT Bounded.")
    
    # TODO WRITE Bound check for Unsafe Set

        
    

