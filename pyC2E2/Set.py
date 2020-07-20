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
    def __init__(self, expressions=None, Eq_Matrix=None):
        """ Construct Set Using either expressions or Eq Matrix(a, b, eq) """
        self.name = "RectangleSet Default Name"
        self.raw_expressions = expressions
        if expressions is not None:
            self.raw_to_matrix()
        elif Eq_Matrix is not None:
            self.aMatrix = Eq_Matrix[0]
            self.bMatrix = Eq_Matrix[1]
            self.eqMatrix = Eq_Matrix[2]
        else:
            self.aMatrix = None
            self.bMatrix = None
            self.eqMatrix = None   

    def __repr__(self):
        return self.raw_expressions 

    def raw_to_matrix(self):
        """ Using SymEq's get_eqn_matrix function """
        # extract variable list (initial set guarantees to include all vars)
        expressions = self.raw_expressions
        # Remove "Mode information"
        pure_expressions = expressions.split(":")[1]
        varList = SymEq.get_var_list(pure_expressions)
        aMatrix, bMatrix, eqMatrix = SymEq.get_eqn_matrix(pure_expressions, varList)
        self.aMatrix = aMatrix
        self.bMatrix = bMatrix
        self.eqMatrix = eqMatrix

    # ERROR CHECKING: Syntax
    def check_initial_set_syntax(self):
        """ Using SymEq's Check Expression Syntax """
        expressions = self.raw_expressions
        if SymEq.check_expression_syntax(expressions, "Initial Set"):
            return True
        else:
            raise Exception("[RectangleSet ERROR]: Syntax NOT valid.")

    # ERROR CHECKING: Boundedness
    def check_initial_set_boundedness(self):
        """ Using SymEq's Check boundedness """
        # Remove "Mode information"
        pure_expressions = self.raw_expressions.split(":")[1]
        if SymEq.check_boundedness(self.aMatrix,
                                    self.bMatrix,
                                    self.eqMatrix,
                                    SymEq.get_var_list(pure_expressions)):
            return True
        else:
            raise Exception("[RectangleSet ERROR]: Set NOT Bounded.")

        
    

