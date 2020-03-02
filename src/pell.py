# -*- coding: utf-8 -*-
"""
Created on Tue Dec 31 10:45:12 2019

@author: Rybakov
"""

import math
import time

#---------------------------------------------------------------------------------------------------
# Matrix.
#---------------------------------------------------------------------------------------------------

class M4:
    """
    Matrix class.
    """

#---------------------------------------------------------------------------------------------------

    def __init__(self, e = [0, 0, 0, 0]):
        """
        Create matrix.

        Arguments:
            e -- elements
        """

        self.E = e

#---------------------------------------------------------------------------------------------------

    def Print(self):
        """
        Print matrix.
        """

        print(str(self.E))

#---------------------------------------------------------------------------------------------------

    def I():
        """
        Create I matrix.

        Result:
            I matrix.
        """

        return M4([1, 0, 0, 1])

#---------------------------------------------------------------------------------------------------

    def L():
        """
        Create left matrix.

        Result:
            Left matrix.
        """

        return M4([1, 0, 1, 1])

#---------------------------------------------------------------------------------------------------

    def R():
        """
        Create right matrix.

        Result:
            Right matrix.
        """

        return M4([1, 1, 0, 1])

#---------------------------------------------------------------------------------------------------

    def Tot(self):
        """
        Return total of matrix - sum of elements.

        Result:
            Sum of elements.
        """

        return sum(self.E)

#---------------------------------------------------------------------------------------------------

    def Mul(self, m4):
        """
        Multiplication:
            self = self * m4.

        Arguments:
            m4 -- matrix.
        """

        m00 = self.E[0] * m4.E[0] + self.E[1] * m4.E[2]
        m01 = self.E[0] * m4.E[1] + self.E[1] * m4.E[3]
        m10 = self.E[2] * m4.E[0] + self.E[3] * m4.E[2]
        m11 = self.E[2] * m4.E[1] + self.E[3] * m4.E[3]

        self.E = [m00, m01, m10, m11]

#---------------------------------------------------------------------------------------------------

    def StepLeft(self):
        """
        Step left.
        """

        a, b, b1, c = self.E[0], self.E[1], self.E[2], self.E[3]

        assert(b == b1)

        self.E = [a + 2 * b + c, b + c, b + c, c]

#---------------------------------------------------------------------------------------------------

    def StepRight(self):
        """
        Step right.
        """

        a, b, b1, c = self.E[0], self.E[1], self.E[2], self.E[3]

        assert(b == b1)

        self.E = [a, a + b, a + b, a + 2 * b + c]

#---------------------------------------------------------------------------------------------------

    def IsEq(self, m4):
        """
        Check if is equal.

        Arguments:
            m4 -- matrix.
        """

        return self.E == m4.E

#---------------------------------------------------------------------------------------------------
# Functions.
#---------------------------------------------------------------------------------------------------

def solve_pell(d):
    """
    Solve pell equation.

    Arguments:
        d -- parameter.
    """

    ini = [1, 0, 0, -d]
    A0, A = M4(ini), M4(ini)
    N, L, R = M4.I(), M4.L(), M4.R()

    while True:

        t = A.Tot()

        if t > 0:
            A.StepLeft()
            N.Mul(L)
        elif t < 0:
            A.StepRight()
            N.Mul(R)
        else:
            raise Exception('t == 0')

        if A.IsEq(A0):
            break

    x, y = N.E[0], N.E[2]

    assert(x * x - d * y * y == 1)

    return (x, y)

#---------------------------------------------------------------------------------------------------

def non_squares(n):
    """
    List of all values non squares.

    Result:
        Non squares.
    """

    li = []

    for i in range(1, n):
        s = math.sqrt(i)
        if s != int(s):
            li.append(i)

    return li

#---------------------------------------------------------------------------------------------------
# Tests.
#---------------------------------------------------------------------------------------------------

if __name__ == '__main__':

    start = time.clock()
    res = [(d, solve_pell(d)) for d in non_squares(100000)]
    stop = time.clock()

    print(res)
    print('Time : ', stop - start)

#---------------------------------------------------------------------------------------------------
