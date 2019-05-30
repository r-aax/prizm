# -*- coding: utf-8 -*-
"""
Oriented edge of functional tree.

Created on Thu May 30 12:28:54 2019

@author: Rybakov
"""

class Edge:
    """
    Oriented edge of functional tree.
    """

#---------------------------------------------------------------------------------------------------
# Constructor.
#---------------------------------------------------------------------------------------------------

    def __init__(self, parent, child):
        """
        Constructor of the edge.

        Arguments:
            parent -- parent,
            child -- child.
            """

        self.Parent = parent
        self.Child = child

        # Properties dictionary.
        self.Dict = {}

#---------------------------------------------------------------------------------------------------
# Properties management.
#---------------------------------------------------------------------------------------------------

    def Set(self, prop, val):
        """
        Set property.

        Arguments:
            prop -- property,
            val -- value.
        """

        self.Dict[prop] = val;

#---------------------------------------------------------------------------------------------------

    def Has(self, prop):
        """
        Check if tree has the property with the given name.

        Arguments:
            prop -- property.

        Result:
            True -- if the tree has the property,
            False -- otherwise.
        """

        return prop in self.Dict

#---------------------------------------------------------------------------------------------------

    def GetWithAlternate(self, prop, alt):
        """
        Get property (in None case get alternate).

        Arguments:
            prop -- property,
            alt -- alternate.

        Result:
            Property value.
        """

        if self.Has(prop):
            return self.Dict[prop]
        else:
            return alt

#---------------------------------------------------------------------------------------------------

    def Get(self, prop):
        """
        Get property.

        Arguments:
            prop -- property.

        Result:
            Property value.
        """

        return self.GetWithAlternate(prop, None)

#---------------------------------------------------------------------------------------------------

    def Is(self, prop, val):
        """
        Check property.

        Arguments:
            prop -- property,
            val -- value.

        Result:
            True -- if property 'prop' is equal to 'val',
            False -- otherwise.
        """

        return self.Get(prop) == val

#---------------------------------------------------------------------------------------------------

    def Del(self, prop):
        """
        Delete property.

        Arguments:
            prop -- property.
        """

        if self.Has(prop):
            self.Dict.pop(prop)

#---------------------------------------------------------------------------------------------------

    def ToString(self):
        """
        Convert to string.

        Result:
            String.
        """

        s = ''

        for p in self.Dict:
            v = self.Get(p)
            
            if p == 'count':
                # Special mean.
                ps = 'x' + str(v)
            else:
                ps = str(p) + ' = ' + str(v)

            if s == '':
                s = ps
            else:
                s = s + ', ' + ps

        return s

#---------------------------------------------------------------------------------------------------
# Tests.
#---------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    e = Edge(1, 2)
    print(e.Source, e.Tree)

#---------------------------------------------------------------------------------------------------
