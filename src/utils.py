# -*- coding: utf-8 -*-
"""
Utils functions.

Created on Mon Jul  2 15:34:59 2018

@author: Rybakov
"""

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
    pass

#---------------------------------------------------------------------------------------------------
