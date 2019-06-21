# -*- coding: utf-8 -*-
"""
Functional.

Created on Wed May 15 10:58:59 2019

@author: Rybakov
"""

import operator

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

def partial_tail2(f, p):
    """
    Partial apply function - set the second argument.

    Arguments:
        f -- function of 2 arguments,
        p -- parameter.

    Result:
        New function of 1 argument.
    """

    return lambda x: f(x, p)

#---------------------------------------------------------------------------------------------------

def partial_tail3(f, p):
    """
    Partial apply function - set the third argument.

    Arguments:
        f -- function of 3 arguments,
        p -- parameter.

    Result:
        New function of 2 arguments.
    """

    return lambda x, y: f(x, y, p)

#---------------------------------------------------------------------------------------------------

def flip(f):
    """
    Flip function.

    Arguments:
        f -- function.

    Result:
        Flipped function.
    """

    return lambda x, y: f(y, x)

#---------------------------------------------------------------------------------------------------

def concat(f, g):
    """
    Concatenate two functions.

    Arguments:
        f -- first function,
        g -- second function.

    Result:
        Concatenated function.
    """

    return lambda x: g(f(x))

#---------------------------------------------------------------------------------------------------

def is_all(l, f):
    """
    Check if all elements of the list satisfy the predicate.

    Arguments:
        l -- list,
        f -- predicate.

    Result:
        True - if all elements of the list satisfy the predicate,
        Fasle - otherwise.
    """

    for e in l:
        if not f(e):
            return False

    return True

#---------------------------------------------------------------------------------------------------

def is_any(l, f):
    """
    Chek if any element of the list satisty the predicate.

    Arguments:
        l -- list,
        f -- predicate.

    Result:
        True - if any element of the list satisfy the predicate,
        False - otherwise.
    """

    return not is_all(l, lambda x: not f(x))

#---------------------------------------------------------------------------------------------------
# Tests.
#---------------------------------------------------------------------------------------------------

if __name__ == '__main__':

    assert partial_tail2(operator.add, 5)(10) == 15
    assert concat(lambda x: x * x, lambda x: x + 5)(2) == 9

#---------------------------------------------------------------------------------------------------
