# -*- coding: utf-8 -*-
"""
Some lists functions.

Created on Thu May 16 14:31:59 2019

@author: Rybakov
"""

#---------------------------------------------------------------------------------------------------
# Functions.
#---------------------------------------------------------------------------------------------------

def descartes_product(a, b):
    """
    Descartes product of two lists.

    Arguments:
        a -- first list,
        b -- second list.

    Result:
        Descartes product.
    """

    lb = len(b)

    return [(a[i // lb], b[i % lb]) for i in range(len(a) * lb)]

#---------------------------------------------------------------------------------------------------

def last(a):
    """
    Last element of the list.

    Arguments:
        a -- list.

    Result:
        Last element.
    """

    if a == []:
        return None
    else:
        return a[len(a) - 1]

#---------------------------------------------------------------------------------------------------

def group(a):
    """
    Group elements of list.

    Arguments:
        a -- list.

    Result:
        New list with groups.
    """

    if a == []:
        return []

    s = sorted(a)
    g = [(s[0], 1)]
    s = s[1:]

    for i in range(len(s)):
        si = s[i]
        (gl, gc) = last(g)
        if gl == si:
            g[len(g) - 1] = (gl, gc + 1)
        else:
            g.append((si, 1))

    return g

#---------------------------------------------------------------------------------------------------

if __name__ == '__main__':

    assert descartes_product(['a', 'b'], [1, 2, 3]) == \
           [('a', 1), ('a', 2), ('a', 3), ('b', 1), ('b', 2), ('b', 3)]
    assert group([1, 3, 2, 1, 3, 2, 3, 1, 3]) == [(1, 3), (2, 2), (3, 4)]
