# -*- coding: utf-8 -*-
"""
Money class.

Created on Tue Feb 19 17:44:14 2019

@author: Rybakov
"""

import utils as ut
import numpy as np
from functools import reduce

class Money:
    """
    Peace of money for financial manipulations.
    """

    HiCost = 100      # Hi part cost (count of "lo" in one "hi")
                      # (1 dollar = 100 cents, etc.).
    DigitsInGroup = 3 # Digits count in group when number represented 
                      # as "1 000 000" (three digits in this case).
    Delim = " "       # Delimiter for parts of numerical value.

#---------------------------------------------------------------------------------------------------
# Constructor and set functions.
#---------------------------------------------------------------------------------------------------

    def __init__(self, v = 0.0):
        """
        Constructor from value (float or string).

        Arguments:
            v -- value.
        """

        # Process value.
        if type(v) is float:
            # Just float value.
            self.SetFloat(v)
        elif type(v) is str:
            # If it is string we must delete all spaces,
            # because we want to process strings like "1 000 000.00".
            self.SetFloat(float(v.replace(Money.Delim, "")))
        else:
            raise ValueError("Wrong money type.")

#---------------------------------------------------------------------------------------------------

    def SetFloat(self, v):
        """
        Set float value.

        Arguments:
            v -- float value.
        """

        self.Amount = int(round(v * float(Money.HiCost)))

#---------------------------------------------------------------------------------------------------

    def FromAmount(amount):
        """
        Create money from amount.

        Arguments:
            amount -- Amount.

        Result:
            Money.
        """

        m = Money()
        m.Amount = amount

        return m

#---------------------------------------------------------------------------------------------------
# Get information.
#---------------------------------------------------------------------------------------------------

    def GetFloat(self):
        """
        Get float value of money.

        Result:
            Value.
        """

        return self.Amount / float(Money.HiCost)

#---------------------------------------------------------------------------------------------------

    def Hi(self):
        """
        High part of money (rubles, dollars, etc.).

        Result:
            High part of money.
        """

        if self.Amount < 0:
            return -((-self.Amount) // Money.HiCost)
        else:
            return self.Amount // Money.HiCost

#---------------------------------------------------------------------------------------------------

    def Lo(self):
        """
        Low part of money (kopecks, cents, etc.).

        Result:
            Low part of money.
        """

        if self.Amount < 0:
            return (-self.Amount) % Money.HiCost
        else:
            return self.Amount % Money.HiCost

#---------------------------------------------------------------------------------------------------

    def HiStr(self):
        """
        High part string representation in form "1 000 000".

        Result:
            High part string representation.
        """

        # Sign.
        hi = self.Hi()
        sign = "-" if hi < 0 else ""

        # Take absulute value.
        if hi < 0:
            hi = -hi

        chopped = ut.chop(str(hi), -3)
        merged = ut.merge(chopped, [Money.Delim] * (len(chopped) - 1))
        return sign + reduce(lambda a, b: a + b, merged)

#---------------------------------------------------------------------------------------------------

    def LoStr(self):
        """
        Low part string representation in form "dd" (two digits).

        Result:
            Low part string representation.
        """

        s = str(self.Lo())

        if len(s) == 1:
            # Add leading zero.
            return "0" + s
        elif len(s) == 2:
            # Pure string.
            return s
        else:
            raise ValueError("Wrong money low value.")

#---------------------------------------------------------------------------------------------------

    def __repr__(self):
        """
        Get money string representation.

        Result:
            String representation.
        """

        # Check for type.
        if not (type(self.Amount) is int):
            raise TypeError("Money amount must be stored as integer.")

        return self.HiStr() + "." + self.LoStr()

#---------------------------------------------------------------------------------------------------

    def Distribute(self, ks):
        """
        Distribute money value between len(ks) money objects according
        with given coefficients.

        Arguments:
            ks -- numpy array of coefficients.

        Result:
            Distributed money (numpy array).
        """

        # Count of coefficients.
        n = len(ks)

        if n == 0:
            # No distribution.
            raise ValueError("No factors for distribute money.")

        if n == 1:
            # Only one factor.
            return self

        # First normalize list.
        nks = ut.npa_norm(ks)

        # Create array for new moneys.
        ms = [0] * n

        # Cycle of initialization array of amounts for new moneys.
        rest = self.Amount
        for i in range(n - 1):
            am = int(round(self.amount * nks[i]))
            rest -= am
            ms[i] = Money.FromAmount(am)

        # The last element calculate from rest.
        ms[n - 1] = Money.FromAmount(rest)

        # Create money objects.
        return ms

#---------------------------------------------------------------------------------------------------

    def Split(self, k):
        """
        Split money.

        Arguments:
            k -- split value (from 0.0 to 1.0).

        Result:
            Tuple of splitted values.
        """

        if not ut.g_in_bounds(k, 0.0, 1.0):
            # Check split factor.
            raise ValueError("Split factor must be in [0.0, 1.0] segment.")

        return self.Distribute(np.array([k, 1.0 - k]))

#---------------------------------------------------------------------------------------------------

    def __add__(self, y):
        """
        Add another money value.

        Arguments:
            y -- added money value.

        Result:
            New money value.
        """

        return Money.FromAmount(self.Amount + y.Amount);

#---------------------------------------------------------------------------------------------------
        
    def __sub__(self, y):
        """
        Sub another money value.

        Arguments:
            y -- subtracted money value.

        Result:
            New money value.
        """

        return Money.FromAmount(self.Amount - y.Amount)

#---------------------------------------------------------------------------------------------------
        
    def __mul__(self, y):
        """
        Multiplication on float value.
        
        Arguments:
            y -- multiplication factor.
        
        Result:
            New money value.
        """
        
        return Money.FromAmount(int(round(self.Amount * y)))

#---------------------------------------------------------------------------------------------------
# Tests.
#---------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    assert True, "stub"

#---------------------------------------------------------------------------------------------------
