# -*- coding: utf-8 -*-
"""
Mathematical functions.

Created on Wed May 15 10:52:29 2019

@author: Rybakov
"""

import fun
import operator
import math

#---------------------------------------------------------------------------------------------------
# Main fucntions.
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
