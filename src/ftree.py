# -*- coding: utf-8 -*-
"""
Functional tree realization.
Tree consists of nodes linked by edges.
Each node is dictionary.

Main data:
    - Type
    - Name
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

    def __init__(self, tp, nm):
        """
        Constructor from type and name.

        Arguments:
            type -- type,
            name -- name.
        """

        self.Dict = {}
        self.Type = tp
        self.Name = nm
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

    def AddChild(self, tp, nm):
        """
        Add new child with given type and name.

        Arguments:
            tp -- type,
            nm -- name.

        Result:
            Added child.
        """

        ch= FTree(tp, nm);
        self.Children.append(ch);
        ch.Parent = self;

        return ch;

#---------------------------------------------------------------------------------------------------
# Print.
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
            own_str = self.Type + " : " + self.Name
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
    earth.AddChild("continent", "Eurasia")
    earth.AddChild("continent", "North America")
    earth.AddChild("continent", "South America")
    earth.AddChild("continent", "Africa")
    earth.AddChild("continent", "Australia")
    earth.AddChild("continent", "Antarctica")

    # Countries.
    with earth.FindElementByTypeName("continent", "Eurasia") as n:
        n.AddChild("country", "Russia")
        n.AddChild("country", "China")
        n.AddChild("country", "Germany")
    with earth.FindElementByTypeName("continent", "North America") as n:
        n.AddChild("country", "USA")
        n.AddChild("country", "Canada")
    with earth.FindElementByTypeName("continent", "South America") as n:
        n.AddChild("country", "Brazil")
        n.AddChild("country", "Argentina")
        n.AddChild("country", "Venezuela")
    with earth.FindElementByTypeName("continent", "Africa") as n:
        n.AddChild("country", "Egypt")
        n.AddChild("country", "RSA")
        n.AddChild("country", "Nigeria")
    with earth.FindElementByTypeName("continent", "Australia") as n:
        n.AddChild("country", "Australia")
    with earth.FindElementByTypeName("continent", "Antarctica") as n:
        # No countries.
        pass

    # Cities.
    with earth.FindElementByTypeName("country", "Russia") as n:
        n.AddChild("city", "Moscow")
        n.AddChild("city", "St. Petersburg")
        n.AddChild("city", "Kazan")
    with earth.FindElementByTypeName("country", "China") as n:
        n.AddChild("city", "Beijing")
        n.AddChild("city", "Shanghai")
    with earth.FindElementByTypeName("country", "Germany") as n:
        n.AddChild("city", "Berlin")
        n.AddChild("city", "Munich")
        n.AddChild("city", "Dresden")
    with earth.FindElementByTypeName("country", "USA") as n:
        n.AddChild("city", "New York")
        n.AddChild("city", "Los Angeles")
        n.AddChild("city", "Chicago")
    with earth.FindElementByTypeName("country", "Canada") as n:
        n.AddChild("city", "Montreal")
    with earth.FindElementByTypeName("country", "Brazil") as n:
        n.AddChild("city", "Rio de Janeiro")
        n.AddChild("city", "San Paulo")
    with earth.FindElementByTypeName("country", "Argentina") as n:
        n.AddChild("city", "Buenos Aires")
    with earth.FindElementByTypeName("country", "Venezuela") as n:
        n.AddChild("city", "Caracas")
    with earth.FindElementByTypeName("country", "Egypt") as n:
        n.AddChild("city", "Cairo")
    with earth.FindElementByTypeName("country", "RSA") as n:
        n.AddChild("city", "Cape Town")
    with earth.FindElementByTypeName("country", "Nigeria") as n:
        n.AddChild("city", "Abuja")
    with earth.FindElementByTypeName("country", "Australia") as n:
        n.AddChild("city", "Sydney")
        n.AddChild("city", "Melbourne")

    earth.PrintTree()

    # Slices.
    for el in earth.SliceChildrenNumbers([1, 0, 2]):
        el.PrintOne()
