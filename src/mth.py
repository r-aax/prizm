# -*- coding: utf-8 -*-
"""
Mathematical functions.

Created on Wed May 15 10:52:29 2019

@author: Rybakov
"""

import fun
import operator

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
# Aggregate functions.
#---------------------------------------------------------------------------------------------------

def avg_arith(vs):
    """
    Average arithmetic.

    Arguments:
        vs -- values.
    """

    return sum(vs) / len(vs)

#---------------------------------------------------------------------------------------------------

def avg_weighted(vs, ws):
    """
    Average value.

    Arguments:
        vs -- values,
        ws -- weights.
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

