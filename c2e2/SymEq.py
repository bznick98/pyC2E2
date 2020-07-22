# SymEq.py Class List:
# SymEq()

import re
from sympy import core as sympy
from scipy import optimize as opt


class SymEq:
    """ Symbolic Equation Library """

    @staticmethod
    def construct_eqn(eqn, is_eq, rationalize):
        """ 
        Construct equation with sympy. Rationalize it if flag is set to True. is_eq is True for DAIs and Actions 
        
        NOTE: sympy.sympify(eqn) converts eqn into a type that can be used inside SymPy
        """

        try:
            if is_eq:
                eqn_split = eqn.split('=')
                lhs, rhs = eqn_split[0], eqn_split[1]
                #print (lhs, rhs)
                eqn = sympy.Eq(sympy.sympify(lhs), sympy.sympify(rhs))
            else:
                eqn = sympy.sympify(eqn)
                if type(eqn) is bool:
                    return eqn
        except:
            raise Exception("[Symbolic Equation ERROR]: Invalid expression.\n")

        if rationalize:
            eqn = SymEq.rationalize(eqn)
        return eqn

    def to_str(self):
        # FIXME self.eqn not exists 
        return str(self.eqn)

    # Given an equation convert all decimals to integers by multiplying by LCM
    @staticmethod
    def rationalize(expr):
        if expr.is_Relational:
            mult = max(SymEq.get_factor(expr.lhs), SymEq.get_factor(expr.rhs))
            return expr.func(mult*expr.lhs, mult*expr.rhs)

    @staticmethod
    def get_factor(expr):
        exp = 0
        if expr.is_Add:
            terms = expr.args
        else:
            terms = [expr]

        for term in terms:
            if term.is_Number:
                exp = max(exp, SymEq.get_term_exp(term))
            for unit in term.args:
                if unit.is_Number:
                    exp = max(exp, SymEq.get_term_exp(unit))

        return 10**exp

    @staticmethod
    def get_term_exp(unit):
        # print("****test****")
        # print(unit)
        text=str(float(unit))
        try:
            idx=text.index('e')
        except ValueError:
            val = len(str(float(unit))) - str((float(unit))).index('.') - 1
            # print(val)
            return val
        exp1 = abs(int(text[idx + 1:len(text)]))
        dec = text[0:idx]
        try:
            exp2 = len(dec)-dec.index('.')-1
        except ValueError:
            # print(exp1)
            return exp1
        val = exp1+exp2
        # print(val)
        return val
        
    @staticmethod
    def get_eqn_matrix(expressions, varList):
        """ Return A and B matrices for expression with x as varList """

        exprs = []
        for expr in expressions.split('&&'):
            if '==' in expr:
                eq = expr.split('==')
                exprs.append(sympy.sympify(eq[0]+'<='+eq[1]))
                exprs.append(sympy.sympify(eq[0]+'>='+eq[1]))
            else:
                sym_expr = sympy.sympify(expr)
                if sym_expr.func==sympy.StrictLessThan:
                    sym_expr = sympy.LessThan(sym_expr.lhs, sym_expr.rhs)
                elif sym_expr.func==sympy.StrictGreaterThan:
                    sym_expr = sympy.GreaterThan(sym_expr.lhs, sym_expr.rhs)
                exprs.append(sym_expr)
        aMatrix=[[0 for v in varList] for expr in exprs]
        bMatrix=[0 for expr in exprs]
        eqMatrix=[]
        for i, expr in enumerate(exprs):
            eqMatrix.append([expr.rel_op])
            expr = expr.lhs-expr.rhs
            expr *= SymEq.get_factor(expr)
            for v in expr.free_symbols:
                aMatrix[i][varList.index(str(v))] = expr.as_coefficients_dict()[v]
                bMatrix[i] = [-expr.as_coefficients_dict()[sympy.numbers.One()]] # FIXME: maybe update sympy version?
        return [aMatrix, bMatrix, eqMatrix]

    # Return whether the initial set is bounded or not
    @staticmethod
    def check_boundedness(a_m, b_m, eq_m, var_list):

        A, B = [], []
        for r, a in enumerate(a_m):
            if eq_m[r][0] == '<=':
                A.append(a)
                B.append(b_m[r][0])
            else:
                A.append([-1 * x for x in a])
                B.append(-1 * b_m[r][0])

        # Solve minimum
        C = [0.0 for v in var_list]
        bounds = [(None, None) for v in var_list]
        for i, v in enumerate(var_list):
            C[i] = 1.0
            min_ = opt.linprog(C, A_ub=A, b_ub=B, bounds=bounds, 
                options=dict(tol=1e-11))  # Added tolerance to eliminate incorrect failures
            if not min_.success:
                return False

            C[i] = 0.0

        return True

    @staticmethod
    def check_expression_syntax(raw_expressions, isInitialSet=False):
        """ Check Initial Set Syntax """
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

        # checking Mode or not checking Mode
        if isInitialSet:
            checking_target = re_is
        else:
            checking_target = re_us

        match = re.match(checking_target, raw_expressions)

        if(match == None):
            return False
        else:
            return True    


    # Represent u**n as u*u*...u for the simulators.
    @staticmethod
    def convert_pow(expr):
        #SymEq.pow_to_mul(expr)
        return str(SymEq.convert_pow_helper(expr))

    @staticmethod
    def convert_pow_helper(expr):
        if not expr.args:
            return expr
        conv_args = [SymEq.convert_pow_helper(arg) for arg in expr.args]
        if expr.is_Pow and expr.exp>0:
            print("this expr is pow:", expr)
            print(expr.base, expr.exp)
            return sympy.Symbol('*'.join(['('+str(conv_args[0])+')']*int(conv_args[1])))
        return expr.func(*conv_args)

    @staticmethod
    def pow_to_mul(expr):
        """
        Convert integer powers in an expression to Muls, like a**2 => a*a.
        """
        pows = list(expr.atoms(sympy.Pow))
        if any(not e.is_Integer for b, e in (i.as_base_exp() for i in pows)):
            raise ValueError("A power contains a non-integer exponent")
        print(pows)
        for b,e in (i.as_base_exp() for i in pows):
            print(b,e)
        repl = zip(pows, (sympy.Mul(*[b]*e,evaluate=False) for b,e in (i.as_base_exp() for i in pows)))
        #print repl
        return expr.subs(repl)

    # Negate guard to construct invariant
    @staticmethod
    def construct_invariant(guard):
        inv_eqn = []
        for eqn in guard.expr:
            print('Eqn is: ' + str(eqn))
            if eqn.func==sympy.LessThan:
                inv = sympy.GreaterThan(*eqn.args)
            elif eqn.func==sympy.GreaterThan:
                inv = sympy.LessThan(*eqn.args)    
            elif eqn.func==sympy.StrictLessThan:
                inv = sympy.StrictGreaterThan(*eqn.args)    
            elif eqn.func==sympy.StrictGreaterThan:
                inv = sympy.StrictLessThan(*eqn.args)
            inv_eqn.append(inv)
        inv_eqn = [str(inv) for inv in inv_eqn]
        return Invariant('||'.join(inv_eqn))

    # Return all free vars in exprs
    @staticmethod
    def vars_used(exprs):
        var_set = set()
        for eqn in exprs:
            var_set = var_set.union(eqn.free_symbols)
        return [str(v) for v in var_set]

    # Return if expr is linear
    @staticmethod
    def is_linear(expr):
        syms = expr.free_symbols
        for x in syms:
            for y in syms:
                try:
                    if not sympy.Eq(sympy.diff(expr, x, y), 0):
                        return False
                except TypeError:
                    return False
        return True

    @staticmethod
    def get_var_list(expressions):
        """ given RAW expressions, extract all variables appeared in the expression """
        varList = []
        for expr in expressions.split('&&'):
            eq = re.split('>|>=|<|<=|==', expr)
            var = eq[0].strip()
            varList.append(var)
        # remove duplicates
        varList = list(set(varList))
        return varList