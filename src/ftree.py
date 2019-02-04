# -*- coding: utf-8 -*-
"""
Functional tree realization.
Tree consists of nodes linked by edges.
Each node is dictionary.

Data is stored in Dict member.
Predefined node keys are:
    - type
    - name

Other data:
    - Children

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
        self.SetType(tp)
        self.SetName(nm)
        self.EmptyChildren()

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
# Basic properties of elements.
#---------------------------------------------------------------------------------------------------

    def GetType(self):
        """
        Get special "type" field value of the element.

        Result:
            Value of the field "type".
        """

        return self.Dict["type"]

#---------------------------------------------------------------------------------------------------

    def SetType(self, tp):
        """
        Set special "type" field of the element.

        Arguments:
            tp -- value.
        """

        self.Dict["type"] = tp

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

        return self.GetType() == tp

#---------------------------------------------------------------------------------------------------

    def GetName(self):
        """
        Get special "name" field value of the element.

        Result:
            Value of the field "name".
        """

        return self.Dict["name"]

#---------------------------------------------------------------------------------------------------

    def SetName(self, nm):
        """
        Set special "name" field of the element.

        Arguments:
            nm -- value.
        """

        self.Dict["name"] = nm

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

        return self.GetName() == nm

#---------------------------------------------------------------------------------------------------

    def GetChildren(self):
        """
        Get special "children" field value of the element.

        Result:
            Value of the field "children".
        """

        return self.Children

#---------------------------------------------------------------------------------------------------
# Elements management.
#---------------------------------------------------------------------------------------------------

    def EmptyChildren(self):
        """
        Empty children list.
        """

        self.Children = []

#---------------------------------------------------------------------------------------------------

    def AddChild(self, child):
        """
        Add new child.

        Arguments:
            child -- new child.

        Result:
            Added child.
        """

        self.GetChildren().append(child)

        return child

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
        print(sh_str + self.GetType() + " : " + self.GetName())

        # Print all children.
        if is_recursive:
            for c in self.GetChildren():
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
            for c in self.GetChildren():
                if c.FindElement(fun) != None:
                    return c

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

        return self.FindElement(lambda x: x.IsName(nm))

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

        return self.FindElement(lambda x: x.IsType(tp) and x.IsName(nm))

#---------------------------------------------------------------------------------------------------
# Slices.
#---------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------
# Tests.
#---------------------------------------------------------------------------------------------------

if __name__ == "__main__":

    print("ftree tests:")

    # Main tree.
    earth = FTree("planet", "Earth")
    
    # Continents.
    earth.AddChild(FTree("continent", "Eurasia"))
    earth.AddChild(FTree("continent", "North America"))
    earth.AddChild(FTree("continent", "South America"))
    earth.AddChild(FTree("continent", "Africa"))
    earth.AddChild(FTree("continent", "Australia"))
    earth.AddChild(FTree("continent", "Antarctica"))

    # Countries.
    with earth.FindElementByTypeName("continent", "Eurasia") as c:
        c.AddChild(FTree("country", "Russia"))
        c.AddChild(FTree("country", "China"))
        c.AddChild(FTree("country", "Germany"))
    with earth.FindElementByTypeName("continent", "North America") as c:
        c.AddChild(FTree("country", "USA"))
        c.AddChild(FTree("country", "Canada"))
    with earth.FindElementByTypeName("continent", "South America") as c:
        c.AddChild(FTree("country", "Brazil"))
        c.AddChild(FTree("country", "Argentina"))
        c.AddChild(FTree("country", "Venezuela"))
    with earth.FindElementByTypeName("continent", "Africa") as c:
        c.AddChild(FTree("country", "Egypt"))
        c.AddChild(FTree("country", "RSA"))
        c.AddChild(FTree("country", "Nigeria"))
    with earth.FindElementByTypeName("continent", "Australia") as c:
        c.AddChild(FTree("country", "Australia"))
    with earth.FindElementByTypeName("continent", "Antarctica") as c:
        # No countries.
        pass

    # Cities.
    with earth.FindElementByTypeName("country", "Russia") as c:
        c.AddChild(FTree("city", "Moscow"))
        c.AddChild(FTree("city", "St. Petersburg"))
        c.AddChild(FTree("city", "Kazan"))
    with earth.FindElementByTypeName("country", "China") as c:
        c.AddChild(FTree("city", "Beijing"))
        c.AddChild(FTree("city", "Shanghai"))
    with earth.FindElementByTypeName("country", "Germany") as c:
        c.AddChild(FTree("city", "Berlin"))
        c.AddChild(FTree("city", "Munich"))
        c.AddChild(FTree("city", "Dresden"))
    with earth.FindElementByTypeName("country", "USA") as c:
        c.AddChild(FTree("city", "New York"))
        c.AddChild(FTree("city", "Los Angeles"))
        c.AddChild(FTree("city", "Chicago"))
    with earth.FindElementByTypeName("country", "Canada") as c:
        c.AddChild(FTree("city", "Montreal"))
    with earth.FindElementByTypeName("country", "Brazil") as c:
        c.AddChild(FTree("city", "Rio de Janeiro"))
        c.AddChild(FTree("city", "San Paulo"))
    with earth.FindElementByTypeName("country", "Argentina") as c:
        c.AddChild(FTree("city", "Buenos Aires"))
    with earth.FindElementByTypeName("country", "Venezuela") as c:
        c.AddChild(FTree("city", "Caracas"))
    with earth.FindElementByTypeName("country", "Egypt") as c:
        c.AddChild(FTree("city", "Cairo"))
    with earth.FindElementByTypeName("country", "RSA") as c:
        c.AddChild(FTree("city", "Cape Town"))
    with earth.FindElementByTypeName("country", "Nigeria") as c:
        c.AddChild(FTree("city", "Abuja"))
    with earth.FindElementByTypeName("country", "Australia") as c:
        c.AddChild(FTree("city", "Sydney"))
        c.AddChild(FTree("city", "Melbourne"))

    earth.PrintTree()
