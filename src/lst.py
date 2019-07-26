# -*- coding: utf-8 -*-
"""
Some lists functions.

Created on Thu May 16 14:31:59 2019

@author: Rybakov
"""

import random
from functools import reduce

#---------------------------------------------------------------------------------------------------
# Functions.
#---------------------------------------------------------------------------------------------------

def rnd_01(n):
    """
    Random list with values from 0 to 1.

    Arguments:
        n -- list length.

    Result:
        Random list.
    """

    return [random.random() for x in range(0, n)]

#---------------------------------------------------------------------------------------------------

def is_flat(a):
    """
    Check if a list is flat.

    Arguments:
        a -- list.

    Result:
        True -- if the list is flat,
        False -- if the list is not flat.
    """

    # Empty list.
    if a == []:
        return True

    # General case.
    return not isinstance(a[0], list) and is_flat(a[1:])

#---------------------------------------------------------------------------------------------------

def flatten(a):
    """
    Flatten list.

    Argument:
        a -- list.

    Result:
        Flattened list.

    Examples:
        flatten([[1]]) -> [1]
        flatten([1, [2, 3]]) -> [1, 2, 3]
    """

    # Empty list.
    if a == []:
        return []

    # Single element.
    if not isinstance(a, list):
        return [a]

    # Recursion.
    return reduce(lambda x, y: x + flatten(y), a, [])

#---------------------------------------------------------------------------------------------------

def chop(s, size = 1):
    """
    Chop list or string.

    Arguments:
        s -- string,
        size -- chop size (if size is positive then the string is chopped
                from the beginning to its end, if size is negative then the
                string is chopped from the end to its beginning).

    Result:
        List of chunks - chopped list.

    Examples:
        chop("123456789", 4) -> ["1234", "5678", "9"]
        chop("123456789", -4) -> ["1", "2345", "6789"]
    """

    # We neeed to rewrite this function because
    # recursion causes core fault.

    c = s
    r = []
    
    if size > 0:
        # Positive chop size - chop from its head.
        while len(c) > size:
            r.append(c[:size])
            c = c[size:]
        if c != []:
            r.append(c)
    elif size < 0:
        # Negative chop size - chop from its end.
        while len(c) > -size:
            r.insert(0, c[size:])
            c = c[:size]
        if c != []:
            r.insert(0, c)            
    else:
        # Chop size must not be zero.
        raise ValueError("Zero chop size.")

    return r

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
        (gl, gc) = g[-1]
        if gl == si:
            g[len(g) - 1] = (gl, gc + 1)
        else:
            g.append((si, 1))

    return g

#---------------------------------------------------------------------------------------------------

def usort(a):
    """
    Sort list and get only unique values.

    Arguments:
        a -- list.

    Result:
        Sorted list with only unique values.
    """

    # Empty list.
    if a == []:
        return []

    # First sort the list.
    s = sorted(a)

    return reduce(lambda r, e: r if r[-1] == e else r + [e], s[1:], [s[0]])

#---------------------------------------------------------------------------------------------------

def mask(li, mask):
    """
    Mask list.

    Arguments:
        li -- list,
        mask -- mask.

    Result:
        Masked list.
    """

    return [e for (e, m) in zip(li, mask) if m == 1]

#---------------------------------------------------------------------------------------------------

def dot(a, b):
    """
    Dot product.

    Arguments:
        a -- first list,
        b -- second list.

    Result:
        Dot product.
    """

    return sum([ai * bi for (ai, bi) in zip(a, b)])

#---------------------------------------------------------------------------------------------------

def merge(a, b):
    """
    Merge two lists.

    Arguments:
        a -- the first list,
        b -- the second list.

    Result:
        Merged list.
    
    Examples:
        merge([1, 1, 1], [2, 2, 2]) -> [1, 2, 1, 2, 1, 2]
        merge([1, 1, 1, 1], [2, 2]) -> [1, 2, 1, 2, 1, 1]
        merge([1, 1], [2, 2, 2, 2]) -> [1, 2, 1, 2, 2, 2]
    """

    # Check the case when one list is empty.
    if a == []:
        return b
    if b == []:
        return a

    # Both lists are non empty.
    return [a[0], b[0]] + merge(a[1 :], b[1 :])

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
# Tests.
#---------------------------------------------------------------------------------------------------

if __name__ == '__main__':

    assert is_flat([1, 2, 3])
    assert not is_flat([1, [2], 3])
    #
    assert flatten(1) == [1]
    assert flatten([[1]]) == [1]
    assert flatten([1, [2, 3]]) == [1, 2, 3]
    #
    assert chop("123456789", 4) == ["1234", "5678", "9"]
    assert chop("123456789", -4) == ["1", "2345", "6789"]
    assert chop([1, 2, 3, 4], 2) == [[1, 2], [3, 4]]
    assert chop([1, 2, 3, 4, 5, 6, 7], 3) == [[1, 2, 3], [4, 5, 6], [7]]
    #
    assert group([1, 3, 2, 1, 3, 2, 3, 1, 3]) == [(1, 3), (2, 2), (3, 4)]
    #
    assert usort([3, 2, 4, 1, 2, 5, 3, 2, 4]) == [1, 2, 3, 4, 5]
    #
    assert mask([1, 2, 3, 4, 5], [0, 0, 0, 0, 0]) == []
    assert mask([1, 2, 3, 4, 5], [1, 1, 0, 0, 1]) == [1, 2, 5]
    #
    assert dot([1, 2, 3], [2, 3, 4]) == 20
    #
    assert merge([1, 1], [2, 2]) == [1, 2, 1, 2]
    assert merge([1], [2, 2, 2]) == [1, 2, 2, 2]
    assert merge([1, 1, 1], [2]) == [1, 2, 1, 1]
    #
    assert descartes_product(['a', 'b'], [1, 2, 3]) == \
           [('a', 1), ('a', 2), ('a', 3), ('b', 1), ('b', 2), ('b', 3)]

#---------------------------------------------------------------------------------------------------
