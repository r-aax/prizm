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

    def Get(self, prop):
        """
        Get property.

        Arguments:
            prop -- property.

        Result:
            Property value.
        """

        if self.Has(prop):
            return self.Dict[prop]
        else:
            return None

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

    def ChildrenNames(self):
        """
        Get all children names.

        Result:
            Children names.
        """

        return [ch.Name for ch in self.Children]

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

    def IsList(self):
        """
        List check.

        Result:
            True -- if is list,
            False -- if is not list.
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

    def AddChildTree(self, t, count = 1):
        """
        Add new child with given type and name.

        Arguments:
            t -- child FTree,
            count -- count (weight).

        Result:
            Added child.
        """

        # Links.
        self.Children.append(t);
        t.Parent = self;

        # Add count.
        if count > 1:
            prop = t.CountStr()
            self.Set(prop, count)

        return t;

#---------------------------------------------------------------------------------------------------

    def AddChildTN(self, tp, nm, descr = "", count = 1):
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

        return self.AddChildTree(t, count);

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

    def CountStr(self):
        """
        Count string for parent.

        Result:
            Count string.
        """

        return self.Type + "_" + self.Name + "_count"

#---------------------------------------------------------------------------------------------------

    def ChildrenCountsStr(self):
        """
        String with children counts.

        Result:
            String with children counts.
        """

        s = ""
        for ch in self.Children:
            prop = ch.CountStr()
            if self.Has(prop):
                s = s + ", " + prop + " = " + str(self.Get(prop))

        # Delete first comma and space.
        if len(s) > 0:
            s = s[2:]

        return s;

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

    def Apply(self, apply_fun, filter_fun = None):
        """
        Apply "apply_fun" function to all elements of the tree.

        Arguments:
            apply_fun -- apply function.
            filter_fun -- additional filter function.
        """

        is_apply = (filter_fun == None) or filter_fun(self)

        # Apply.
        if is_apply:
            apply_fun(self)
            
        # Children.
        for ch in self.Children:
            ch.Apply(apply_fun, filter_fun)

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
