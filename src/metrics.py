# -*- coding: utf-8 -*-
"""
Metrics.

Created on Sun May 19 23:06:35 2019

@author: Rybakov
"""

import math
import mth

#---------------------------------------------------------------------------------------------------

def lp_norm(a, b, p):
    """
    l_p norm.

    Arguments:
        a -- first point,
        b -- second point.

    Result:
        Norm result.
    """

    return pow(sum([pow(abs(ai - bi), p) for (ai, bi) in zip(a, b)]), 1.0 / p)

#---------------------------------------------------------------------------------------------------

def sup_norm(a, b):
    """
    Supremum norm.

    Arguments:
        a -- first point,
        b -- second point.

    Result:
        Norm result.
    """

    return max([abs(ai - bi) for (ai, bi) in zip(a, b)])

#---------------------------------------------------------------------------------------------------

def jeffreys_matsushita(a, b):
    """
    Jeffreys-Matsushita.

    Aarguments:
        a -- first point,
        b -- second point.

    Result:
        Near function.
    """

    sqrt_red = lambda x: mth.sign(x) * math.sqrt(abs(x))
    sqrt_diff = lambda x, y: sqrt_red(x) - sqrt_red(y)

    return math.sqrt(sum([pow(sqrt_diff(ai, bi), 2.0) for (ai, bi) in zip(a, b)]))

#---------------------------------------------------------------------------------------------------

def div_coef(a, b):
    """
    Divergence coefficient.

    Arguments:
        a -- first point,
        b -- secod point.

    Result:
        Coefficient value.
    """

    k = 1.0 / len(a)
    f = lambda ai, bi: (ai - bi) / (abs(ai) + abs(bi))

    return math.sqrt(k * sum([pow(f(ai, bi), 2.0) for (ai, bi) in zip(a, b)]))

#---------------------------------------------------------------------------------------------------