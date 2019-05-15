# -*- coding: utf-8 -*-
"""
Functional.

Created on Wed May 15 10:58:59 2019

@author: Rybakov
"""

#---------------------------------------------------------------------------------------------------

def zipwith(a, b, f):
    """
    Zip two list with a function.

    Arguments:
        a -- first list,
        b -- seconds list,
        f -- function.

    Result:
        Zipped list.
    """

    return [f(ai, bi) for (ai, bi) in zip(a, b)]

#---------------------------------------------------------------------------------------------------

def unzip(ab):
    """
    Unzip list of tuples to tuple of lists.

    Arguments:
        ab -- list of tuples.

    Result:
        Tuple of lists.
    """

    return ([ai for (ai, bi) in ab], [bi for (ai, bi) in ab])

#---------------------------------------------------------------------------------------------------
