# -*- coding: utf-8 -*-
"""
Mathematical functions.

Created on Wed May 15 10:52:29 2019

@author: Rybakov
"""

import fun
import operator

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
