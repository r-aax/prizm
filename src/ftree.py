# -*- coding: utf-8 -*-
"""
Functional tree realization.
Tree consists of nodes linked by edges.
Each node is dictionary.

Main data:
    - Type
    - Name
    - Descr
    - Children
    - Parent

Other data is stored in Dict member.

Created on Tue Jan 15 16:05:25 2019

@author: Rybakov
"""

class FTree:
    """
    Functional tree.
    """

#---------------------------------------------------------------------------------------------------
# Constructor.
#---------------------------------------------------------------------------------------------------

    def __init__(self, tp, nm, descr = ""):
        """
        Constructor from type and name.

        Arguments:
            type -- type,
            name -- name,
            descr -- description.
        """

        self.Dict = {}
        self.Type = tp
        self.Name = nm
        self.Descr = descr
        self.Children = []
        self.Parent = None

#---------------------------------------------------------------------------------------------------
# Maintenance.
#---------------------------------------------------------------------------------------------------

    def __enter__(self):
        """
        Function for "with ... as" context.
        """

        return self

#---------------------------------------------------------------------------------------------------

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Method for "with ... as" context.

        Arguments:
            exc_type -- exception type,
            exc_val -- exception value,
            exc_tb -- exception traceback.
        """

        pass

#---------------------------------------------------------------------------------------------------
# Manipulation with basic fields of elements.
#---------------------------------------------------------------------------------------------------

    def IsType(self, tp):
        """
        Check if the element has given type.

        Arguments:
            tp -- value of the field "type".

        Result:
            True -- if the element has given type,
            False -- otherwise.
        """

        return self.Type == tp

#---------------------------------------------------------------------------------------------------

    def IsName(self, nm):
        """
        Check if the element has given name.

        Arguments:
            nm -- value of the field "name".

        Result:
            True -- if the element has given name,
            False -- otherwise.
        """

        return self.Name == nm

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

    def SetOuter(self, ch, prop, val):
        """
        Set outer property for child.

        Arguments:
            ch -- child,
            prop -- property name,
            val -- value.
        """

        self.Set(ch.OuterPropertyStr(prop), val);

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

    def HasOuter(self, ch, prop):
        """
        Check if tree has the outer property with the given name.

        Arguments:
            ch -- child,
            prop -- property.

        Result:
            True -- if the tree has the property,
            False -- otherwise.
        """

        return self.Has(ch.OuterPropertyStr(prop))

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

    def GetOuterWithAlternate(self, ch, prop, alt):
        """
        Get outer property (in None case get alternate).

        Arguments:
            ch -- child,
            prop -- property,
            alt -- alternate.

        Result:
            Property value.
        """

        return self.GetWithAlternate(ch.OuterPropertyStr(prop), alt)

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

    def GetOuter(self, ch, prop):
        """
        Get outer property.

        Arguments:
            ch -- child,
            prop -- property.

        Result:
            Property value.
        """

        return self.Get(ch.OuterPropertyStr(prop))

#---------------------------------------------------------------------------------------------------

    def Is(self, prop, val):
        """
        Check property.

        Arguments:
            prop -- property,
            val -- value.

        Result:
            True -- if property "prop" is equal to "val",
            False -- otherwise.
        """

        return self.Get(prop) == val

#---------------------------------------------------------------------------------------------------

    def IsOuter(self, ch, prop, val):
        """
        Check outer property.

        Arguments:
            ch -- child,
            prop -- property,
            val -- value.

        Result:
            True -- if property "prop" is equal to "val",
            False -- otherwise.
        """

        return self.Is(ch.OuterPropertyStr(prop), val)

#---------------------------------------------------------------------------------------------------
# Properties.
#---------------------------------------------------------------------------------------------------

    def IsRoot(self):
        """
        Root check.

        Result:
            True -- if is root,
            False -- if is not root.
        """

        return self.Parent == None

#---------------------------------------------------------------------------------------------------

    def IsLeaf(self):
        """
        List check.

        Result:
            True -- if is leaf,
            False -- if is not leaf.
        """

        return self.Children == []

#---------------------------------------------------------------------------------------------------

    def ChildrenCount(self):
        """
        Count of children.

        Result:
            Count of children.
        """

        return self.Children.count;

#---------------------------------------------------------------------------------------------------

    def Level(self):
        """
        Node level.

        Result:
            Level.
        """

        if self.IsRoot():
            return 0
        else:
            return 1 + self.Parent.Level()

#---------------------------------------------------------------------------------------------------
# Elements management.
#---------------------------------------------------------------------------------------------------

    def AddChildTree(self, t):
        """
        Add new child with given type and name.

        Arguments:
            t -- child FTree.

        Result:
            Added child.
        """

        # Links.
        self.Children.append(t);
        t.Parent = self;

        return t;

#---------------------------------------------------------------------------------------------------

    def AddChildTN(self, tp, nm, descr = ""):
        """
        Add new child with given type and name.

        Arguments:
            tp -- type,
            nm -- name,
            descr -- description.

        Result:
            Added child.
        """

        t = FTree(tp, nm, descr);

        return self.AddChildTree(t);

#---------------------------------------------------------------------------------------------------
# Print.
#---------------------------------------------------------------------------------------------------

    def BaseStr(self):
        """
        Basic string.

        Result:
            Basic string.
        """

        s = "[" + self.Type + ":" + self.Name + "]"

        if self.Descr != "":
            s = s +" " + self.Descr

        return s

#---------------------------------------------------------------------------------------------------

    def OuterPropertyStr(self, prop):
        """
        Generate outer property string.

        Arguments:
            prop -- property.

        Result:
            Ouuter property string.
        """

        return self.Type + "_" + self.Name + "_" + prop

#---------------------------------------------------------------------------------------------------

    def PropertiesStr(self):
        """
        String with all properties.

        Result:
            Properties string.
        """

        s = ""

        # Add all properties.
        for p in self.Dict:
            v = self.Get(p)
            if not (type(v) is type(lambda x: x)):
                s = s + ", " + p + " = " + str(self.Get(p))

        # Delete first 2 symbols.
        if len(s) > 2:
            s = s[2:]

        return s

#---------------------------------------------------------------------------------------------------

    def Print(self, level, is_recursive):
        """
        Print.

        Arguments:
            level -- tree level,
            is_recursive -- recursive print is needed or not.
        """

        if level == 0:

            # Print "#" ony in recursive mode.
            if is_recursive:
                sh_str = "# " # root
            else:
                sh_str = ""

        else:
            sh = 4 * (level - 1)
            sh_str = (" " * sh) + "|--->"

        # Print current element.
        to_string_fun = self.Get("to_string_fun")
        if to_string_fun != None:
            own_str = to_string_fun(self)
        else:
            own_str = self.BaseStr()
        print(sh_str + own_str)

        # Print all children.
        if is_recursive:
            for c in self.Children:
                c.Print(level + 1, True)

#---------------------------------------------------------------------------------------------------

    def PrintOne(self):
        """
        Print single element.
        """

        self.Print(0, False)

#---------------------------------------------------------------------------------------------------

    def PrintTree(self):
        """
        Recursive print of tree.
        """

        self.Print(0, True)

#---------------------------------------------------------------------------------------------------
# Find elements.
#---------------------------------------------------------------------------------------------------

    def FindElement(self, fun):
        """
        Find element in the tree.

        Arguments:
            fun -- search function.

        Result:
            Found element, if it is in th tree,
            None, if element is not found.
        """

        if fun(self):
            return self
        else:

            # Check all children.
            for c in self.Children:
                f = c.FindElement(fun)
                if f != None:
                    return f

            # Nothing is found.
            return None

#---------------------------------------------------------------------------------------------------

    def FindElementByName(self, nm):
        """
        Find element by its name.

        Arguments:
            nm -- name.

        Result:
            Found element, if it is in the tree,
            None, if element is not found.
        """

        return self.FindElement(lambda t: t.IsName(nm))

#---------------------------------------------------------------------------------------------------

    def FindElementByTypeName(self, tp, nm):
        """
        Find element by its type and name.

        Arguments:
            tp -- type,
            nm -- name.

        Result:
            Found element, if it is in the tree,
            None, if element is not found.
        """

        return self.FindElement(lambda t: t.IsType(tp) and t.IsName(nm))

#---------------------------------------------------------------------------------------------------
# Main functional actions.
#---------------------------------------------------------------------------------------------------

    def Apply(self, apply_fun, filter_fun = None, is_downward = True):
        """
        Apply "apply_fun" function to all elements of the tree.

        Arguments:
            apply_fun -- apply function,
            filter_fun -- additional filter function,
            is_downward -- from up to down.
        """

        is_apply = (filter_fun == None) or filter_fun(self)

        if is_downward:
            # Apply.
            if is_apply:
                apply_fun(self)

            # Children.
            for ch in self.Children:
                ch.Apply(apply_fun, filter_fun, is_downward)
        else:
            # Children.
            for ch in self.Children:
                ch.Apply(apply_fun, filter_fun, is_downward)

            # Apply.
            if is_apply:
                apply_fun(self)

#---------------------------------------------------------------------------------------------------

    def ApplyUpward(self, apply_fun, filter_fun = None):
        """
        Apply "apply_fun" function to all elements of the tree.

        Arguments:
            apply_fun -- apply function,
            filter_fun -- additional filter function.
        """

        self.Apply(apply_fun, filter_fun, False)

#---------------------------------------------------------------------------------------------------

    def ApplyDownward(self, apply_fun, filter_fun = None):
        """
        Apply "apply_fun" function to all elements of the tree.

        Arguments:
            apply_fun -- apply function,
            filter_fun -- additional filter function.
        """

        self.Apply(apply_fun, filter_fun, True)

#---------------------------------------------------------------------------------------------------

    def FoldDepth(self, fun, acc):
        """
        Fold tree while depth-first search.

        Arguments:
            fun -- function,
            acc -- accumulator.

        Result:
            Accumulator value after fold.
        """

        v = fun(self, acc)

        # Now add to accumulator all children results.
        for c in self.Children:
            v = c.FoldDepth(fun, v)

        return v

#---------------------------------------------------------------------------------------------------
# Slices.
#---------------------------------------------------------------------------------------------------

    def SliceLevel(self, level):
        """
        Get level slice.

        Arguments:
            level -- number of level.

        Result:
            Slice.
        """

        return self.FoldDepth(lambda t, a: a + [t] if (t.Level() == level) else a, [])

#---------------------------------------------------------------------------------------------------

    def SliceChildNum(self, n):
        """
        Slice on child number.

        Arguments:
            n -- child number.

        Result:
            Slice.
        """

        s = [self]
        children_len = len(self.Children)

        if (n >= -children_len) and (n < children_len):
            s = s + self.Children[n].SliceChildNum(n)

        return s

#---------------------------------------------------------------------------------------------------

    def SliceLeft(self):
        """
        Slice left path.

        Return:
            Slice.
        """

        return self.SliceChildNum(0)

#---------------------------------------------------------------------------------------------------

    def SliceRight(self):
        """
        Slice left path.

        Return:
            Slice.
        """

        return self.SliceChildNum(-1)

#---------------------------------------------------------------------------------------------------

    def SliceChildrenNumbers(self, ns):
        """
        Slice on children numbers.

        Arguments:
            ns -- children numbers.

        Result:
            Slice.
        """

        s = [self]
        children_len = len(self.Children)

        if (ns != []):
            h = ns[0]
            if (h >= -children_len) and (h < children_len):
                s = s + self.Children[h].SliceChildrenNumbers(ns[1:])

        return s

#---------------------------------------------------------------------------------------------------
# Gather tactics.
#---------------------------------------------------------------------------------------------------

    def GatherTacticSumWithCount(self, prop):
        """
        Gather all 'prop' values from children, multiply on 'count' and
        set to 'prop' of parent.

        Arguments:
            prop -- property.
        """

        if not self.Has(prop):
            self.Set(prop,
                     sum([ch.Get(prop) * self.GetOuterWithAlternate(ch, "count", 1.0)
                          for ch in self.Children]))

#---------------------------------------------------------------------------------------------------
# Tests.
#---------------------------------------------------------------------------------------------------

if __name__ == "__main__":

    print("ftree tests:")

    # Main tree.
    earth = FTree("planet", "Earth")
    
    # Continents.
    earth.AddChildTN("continent", "Eurasia")
    earth.AddChildTN("continent", "North America")
    earth.AddChildTN("continent", "South America")
    earth.AddChildTN("continent", "Africa")
    earth.AddChildTN("continent", "Australia")
    earth.AddChildTN("continent", "Antarctica")

    # Countries.
    with earth.FindElementByTypeName("continent", "Eurasia") as n:
        n.AddChildTN("country", "Russia")
        n.AddChildTN("country", "China")
        n.AddChildTN("country", "Germany")
    with earth.FindElementByTypeName("continent", "North America") as n:
        n.AddChildTN("country", "USA")
        n.AddChildTN("country", "Canada")
    with earth.FindElementByTypeName("continent", "South America") as n:
        n.AddChildTN("country", "Brazil")
        n.AddChildTN("country", "Argentina")
        n.AddChildTN("country", "Venezuela")
    with earth.FindElementByTypeName("continent", "Africa") as n:
        n.AddChildTN("country", "Egypt")
        n.AddChildTN("country", "RSA")
        n.AddChildTN("country", "Nigeria")
    with earth.FindElementByTypeName("continent", "Australia") as n:
        n.AddChildTN("country", "Australia")
    with earth.FindElementByTypeName("continent", "Antarctica") as n:
        # No countries.
        pass

    # Cities.
    with earth.FindElementByTypeName("country", "Russia") as n:
        n.AddChildTN("city", "Moscow")
        n.AddChildTN("city", "St. Petersburg")
        n.AddChildTN("city", "Kazan")
    with earth.FindElementByTypeName("country", "China") as n:
        n.AddChildTN("city", "Beijing")
        n.AddChildTN("city", "Shanghai")
    with earth.FindElementByTypeName("country", "Germany") as n:
        n.AddChildTN("city", "Berlin")
        n.AddChildTN("city", "Munich")
        n.AddChildTN("city", "Dresden")
    with earth.FindElementByTypeName("country", "USA") as n:
        n.AddChildTN("city", "New York")
        n.AddChildTN("city", "Los Angeles")
        n.AddChildTN("city", "Chicago")
    with earth.FindElementByTypeName("country", "Canada") as n:
        n.AddChildTN("city", "Montreal")
    with earth.FindElementByTypeName("country", "Brazil") as n:
        n.AddChildTN("city", "Rio de Janeiro")
        n.AddChildTN("city", "San Paulo")
    with earth.FindElementByTypeName("country", "Argentina") as n:
        n.AddChildTN("city", "Buenos Aires")
    with earth.FindElementByTypeName("country", "Venezuela") as n:
        n.AddChildTN("city", "Caracas")
    with earth.FindElementByTypeName("country", "Egypt") as n:
        n.AddChildTN("city", "Cairo")
    with earth.FindElementByTypeName("country", "RSA") as n:
        n.AddChildTN("city", "Cape Town")
    with earth.FindElementByTypeName("country", "Nigeria") as n:
        n.AddChildTN("city", "Abuja")
    with earth.FindElementByTypeName("country", "Australia") as n:
        n.AddChildTN("city", "Sydney")
        n.AddChildTN("city", "Melbourne")

    earth.PrintTree()

    # Slices.
    for el in earth.SliceChildrenNumbers([1, 0, 2]):
        el.PrintOne()
