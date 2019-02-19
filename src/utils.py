# -*- coding: utf-8 -*-
"""
Utils functions.

Created on Mon Jul  2 15:34:59 2018

@author: Rybakov
"""

import random
from functools import reduce

#---------------------------------------------------------------------------------------------------
# Constants.
#---------------------------------------------------------------------------------------------------

class Consts:
    """
    Class for constants.
    """

    # Basic constants.
    TEN      = 10   # Ten constant.
    HUNDRED  = 100  # Hundred constant.
    THOUSAND = 1000 # Thousand constant.

    # Time constants.
    YEAR_MONTHS    =  12                   # Months count in a year.
    YEAR_DAYS      = 365                   # Days count in a not-leap year.
    WEEK_DAYS      =   7                   # Days count in a week.
    DAY_HOURS      =  24                   # Hours count in a day.
    YEAR_HOURS     = YEAR_DAYS * DAY_HOURS # Hours count in a year
    HOUR_MINUTES   =  60                   # Minutes count in a hour.
    MINUTE_SECONDS =  60                   # Seconds count in a minute.

#---------------------------------------------------------------------------------------------------
# Counter.
#---------------------------------------------------------------------------------------------------

class Cnt:
    """
    Class for counter.
    """

#---------------------------------------------------------------------------------------------------

    def __init__(self):
        """
        Constructor.
        """

        self.v = 0

#---------------------------------------------------------------------------------------------------

    def Inc(self, inc = 1):
        """
        Increment value.

        Arguments:
            inc -- increment value.
        """

        self.v += inc

#---------------------------------------------------------------------------------------------------
# General functions.
#---------------------------------------------------------------------------------------------------

def g_in_bounds(x, lo, hi):
    """
    Check if value is placed in bounds.

    Arguments:
        x -- value,
        lo -- low bound,
        hi -- high bound.

    Result:
        True -- is value is in bounds,
        False -- otherwise.
    """

    return (x >= lo) and (x <= hi)

#---------------------------------------------------------------------------------------------------
# Lists (names start with li_*).
#---------------------------------------------------------------------------------------------------

def li_merge(a, b):
    """
    Merge two lists.

    Arguments:
        a -- the first list,
        b -- the second list.

    Result:
        Merged list.
    
    Examples:
        li_merge([1, 1, 1], [2, 2, 2]) -> [1, 2, 1, 2, 1, 2]
        li_merge([1, 1, 1, 1], [2, 2]) -> [1, 2, 1, 2, 1, 1]
        li_merge([1, 1], [2, 2, 2, 2]) -> [1, 2, 1, 2, 2, 2]
    """

    # Check the case when one list is empty.
    if a == []:
        return b
    if b == []:
        return a

    # Both lists are non empty.
    return [a[0], b[0]] + li_merge(a[1 :], b[1 :])

#---------------------------------------------------------------------------------------------------

def li_flatten(a):
    """
    Flatten list.

    Argument:
        a -- list.

    Result:
        Flattened list.

    Examples:
        li_flatten([[1]]) -> [1]
        li_flatten([1, [2, 3]]) -> [1, 2, 3]
    """

    # Empty list.
    if a == []:
        return []

    # Single element.
    if not isinstance(a, list):
        return [a]

    # Recursion.
    return reduce(lambda x, y: x + li_flatten(y), a, [])

#---------------------------------------------------------------------------------------------------

def li_is_flat(a):
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
    return not isinstance(a[0], list) and li_is_flat(a[1:])

#---------------------------------------------------------------------------------------------------

def li_last(a):
    """
    Get last element of the list.

    Arguments:
        a -- list.

    Result:
        The last element.
    """

    return None if a == [] else a[len(a) - 1]

#---------------------------------------------------------------------------------------------------

def li_sort_uniq(a):
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

    return reduce(lambda r, e: r if li_last(r) == e else r + [e], s[1:], [s[0]])

#---------------------------------------------------------------------------------------------------

def li_rnd_01(n):
    """
    Random list with values from 0 to 1.

    Arguments:
        n -- list length.

    Result:
        Random list.
    """

    return [random.random() for x in range(0, n)]

#---------------------------------------------------------------------------------------------------

def li_chop(s, size = 1):
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
        li_chop("123456789", 4) -> ["1234", "5678", "9"]
        li_chop("123456789", -4) -> ["1", "2345", "6789"]
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
# Numpy arrays (names start with npa_*).
#---------------------------------------------------------------------------------------------------

def npa_norm(a):
    """
    Normalize numerical numpy array.

    Arguments:
        a -- array.

    Result:
        Normalized array.
    """

    return a / sum(a)

#---------------------------------------------------------------------------------------------------
# Tests.
#---------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    assert li_merge([1, 1], [2, 2]) == [1, 2, 1, 2], "li_merge fault 01"
    assert li_merge([1], [2, 2, 2]) == [1, 2, 2, 2], "li_merge fault 02"
    assert li_merge([1, 1, 1], [2]) == [1, 2, 1, 1], "li_merge fault 03"
    #
    assert li_flatten([[1]]) == [1], "li_flatten fault 01"
    assert li_flatten([1, [2, 3]]) == [1, 2, 3], "li__flatten fault 02"
    #
    assert li_chop("123456789", 4) == ["1234", "5678", "9"], "li_chop fault 01"
    assert li_chop("123456789", -4) == ["1", "2345", "6789"], "li_chop fault 02"

#---------------------------------------------------------------------------------------------------
