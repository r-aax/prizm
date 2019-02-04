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

        print("XXX")

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
    jscc = FTree("resources", "JSCC")

    # Adding supercomputers.
    jscc.AddChild(FTree("supercomputer", "100K"))
    jscc.AddChild(FTree("supercomputer", "10P"))

    # 100K segments.
    with jscc.FindElementByTypeName("supercomputer", "100K") as k100:
        k100.AddChild(FTree("segment", "100K"))

    # 10P segments.
    with jscc.FindElementByName("10P") as p10:
        p10.AddChild(FTree("segment", "Tornado"))
        p10.AddChild(FTree("segment", "Petastream"))
        p10.AddChild(FTree("segment", "Haswell"))
        p10.AddChild(FTree("segment", "Broadwell"))
        p10.AddChild(FTree("segment", "KNL"))
        p10.AddChild(FTree("segment", "Skylake"))

    jscc.PrintTree()
