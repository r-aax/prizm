# -*- coding: utf-8 -*-
"""
Mathematical functions.

Created on Wed May 15 10:52:29 2019

@author: Rybakov
"""

import fun
import operator
import math
import lst

#---------------------------------------------------------------------------------------------------
# Matrix class.
#---------------------------------------------------------------------------------------------------

class Matrix:

#---------------------------------------------------------------------------------------------------

    def __init__(self, rows, cols):
        """
        Constructor.

        Arguments:
            rows -- rows count,
            cols -- columns count.
        """

        self.Rows = rows
        self.Cols = cols

        self.Items = [0.0] * rows
        for i in range(rows):
            self.Items[i] = [0.0] * cols

#---------------------------------------------------------------------------------------------------

    def __getitem__(self, i):
        """
        Redefine [].

        Arguments:
            i -- index.

        Result:
            Object.
        """

        return self.Items[i]

#---------------------------------------------------------------------------------------------------

    def __mul__(self, v):
        """
        Multiply on vector.

        Arguments:
            v -- vector.

        Result:
            Result.
        """

        assert isinstance(v, list)
        assert len(v) == self.Cols

        return [lst.dot(row, v) for row in self.Items]

#---------------------------------------------------------------------------------------------------

    def Col(self, j):
        """
        Get column.

        Arguments:
            j -- column index.

        Result:
            Column.
        """

        return [r[j] for r in self.Items]

#---------------------------------------------------------------------------------------------------

    def __rmul__(self, v):
        """
        Multiply on vector (right side).

        Arguments:
            v -- vector.

        Result:
            Result.
        """

        return [lst.dot(v, self.Col(j)) for j in range(self.Cols)]

#---------------------------------------------------------------------------------------------------

    def Set(self, li):
        """
        Set matrix elements.

        Arguments:
            li -- list.
        """

        row, col = 0, 0
        for el in li:
            self[row][col] = el
            col += 1
            if col == self.Cols:
                row += 1
                col = 0

#---------------------------------------------------------------------------------------------------

    def Elements(self):
        """
        Elements count.

        Result:
            Elements count.
        """

        return self.Rows * self.Cols

#---------------------------------------------------------------------------------------------------

    def E(n):
        """
        E matrix.

        Arguments:
            n -- matrix size.

        Result:
            Matrix.
        """

        m = Matrix(n, n)

        for i in range(n):
            m[i][i] = 1.0

        return m

#---------------------------------------------------------------------------------------------------

    def Copy(self):
        """
        Get copy matrix.

        Result:
            Copy.
        """

        rows, cols = self.Rows, self.Cols
        m = Matrix(rows, cols)

        for i in range(rows):
            for j in range(cols):
                m[i][j] = self[i][j]

        return m

#---------------------------------------------------------------------------------------------------

    def IsQuad(self):
        """
        Check if is quad mastrix.

        Result:
            True -- if is quad,
            False -- otherwise.
        """

        return self.Rows == self.Cols

#---------------------------------------------------------------------------------------------------

    def SwapRows(self, i, j):
        """
        Swap rows.

        Arguments:
            i -- first index,
            j -- second index.
        """

        if i == j:
            return

        self.Items[i], self.Items[j] = self.Items[j], self.Items[i]

#---------------------------------------------------------------------------------------------------

    def ScaleRow(self, i, k):
        """
        Scale row.

        Arguments:
            i -- row index,
            k -- coefficient.
        """

        for c in range(self.Cols):
            self[i][c] *= k

#---------------------------------------------------------------------------------------------------

    def AddMulRow(self, i, j, k):
        """
        Add to j-th row i-th row multiplied on k.

        Arguments:
            i -- first row index,
            j -- second row index,
            k -- coefficient.
        """

        for c in range(self.Cols):
            self[j][c] += self[i][c] * k

#---------------------------------------------------------------------------------------------------

    def Inverted(self):
        """
        Get inverted matrix.

        Result:
            Inverted matrix.
        """

        assert self.IsQuad()

        n = self.Rows
        a = self.Copy()
        b = Matrix.E(n)

        for i in range(n):

            # i-th row

            # Lead row.
            max_j, max_v = i, abs(a[i][i])
            for j in range(i + 1, n):
                if abs(a[j][i]) > max_v:
                    max_j, max_v = j, abs(a[j][i])
            a.SwapRows(i, max_j)
            b.SwapRows(i, max_j)

            # Scale row.
            k = 1.0 / a[i][i]
            a.ScaleRow(i, k)
            b.ScaleRow(i, k)

            # Add/Mul.
            for j in range(n):
                if (j != i):
                    k = -a[j][i]
                    a.AddMulRow(i, j, k)
                    b.AddMulRow(i, j, k)

        return b

#---------------------------------------------------------------------------------------------------
# Other fucntions.
#---------------------------------------------------------------------------------------------------

def frac(v):
    """
    Fractional part of float.

    Arguments:
        v -- value.

    Result:
        Fractional part.
    """

    return v - int(v)

#---------------------------------------------------------------------------------------------------

def sign(v):
    """
    Sign ov value.

    Arguments:
        v -- value.

    Result:
        Sign.
    """

    if v > 0:
        return 1.0
    elif v < 0:
        return -1.0
    else:
        return 0.0

#---------------------------------------------------------------------------------------------------

def hypot(x, y, z = 0.0):
    """
    Hypotenuze.

    Arguments:
        x -- X size,
        y -- Y size,
        z -- Z size.

    Result:
        Hypot value.
    """

    return math.sqrt(x * x + y * y + z * z)

#---------------------------------------------------------------------------------------------------
# Aggregate functions.
#---------------------------------------------------------------------------------------------------

def avg_arith(vs):
    """
    Average arithmetic.

    Arguments:
        vs -- values.

    Result:
        Arithmetic average value.
    """

    return sum(vs) / len(vs)

#---------------------------------------------------------------------------------------------------

def avg_weighted(vs, ws):
    """
    Average value.

    Arguments:
        vs -- values,
        ws -- weights.

    Result:
        Weighted average value.
    """

    return sum(fun.zipwith(vs, ws, operator.mul)) / sum(ws)

#---------------------------------------------------------------------------------------------------

def min_with_index(a):
    """
    Minimum value and its index.

    Arguments:
        a -- list.

    Result:
        Minimum value and index.
    """

    if a == []:
        return (None, -1)

    cm, ci = a[0], 0
    for i in range(1, len(a)):
        if a[i] < cm:
            cm, ci = a[i], i

    return (cm, ci)

#---------------------------------------------------------------------------------------------------

def max_with_index(a):
    """
    Maximum value and its index.

    Arguments:
        a -- list.

    Result:
        Maximum value and index.
    """

    if a == []:
        return (None, -1)

    cm, ci = a[0], 0
    for i in range(1, len(a)):
        if a[i] > cm:
            cm, ci = a[i], i

    return (cm, ci)

#---------------------------------------------------------------------------------------------------

def sigmoid(x):
    """
    Sigmoid function.

    Arguments:
        x -- argument.

    Result:
        Result.
    """

    return 1.0 / (1.0 + math.exp(-x))

#---------------------------------------------------------------------------------------------------

def sigmoid_der(x):
    """
    Derivative of sigmoid function.

    Arguments:
        x -- argument.

    Result:
        Result.
    """

    return sigmoid(x) * (1.0 - sigmoid(x))

#---------------------------------------------------------------------------------------------------
# Tests.
#---------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    pass

#---------------------------------------------------------------------------------------------------
