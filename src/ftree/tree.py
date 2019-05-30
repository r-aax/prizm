# -*- coding: utf-8 -*-
"""
Functional tree.

Created on Thu May 30 12:28:38 2019

@author: Rybakov
"""

from ftree.edge import Edge

class Tree:
    """
    Functional tree.
    """

#---------------------------------------------------------------------------------------------------
# Constructor.
#---------------------------------------------------------------------------------------------------

    def __init__(self, tp, nm, descr = ''):
        """
        Constructor from type and name.

        Arguments:
            type -- type,
            name -- name,
            descr -- description.
        """

        self.Type = tp
        self.Name = nm
        self.Descr = descr

        # Properties dictionary.
        self.Dict = {}

        # Array of children edges.
        self.Succs = []

        # Parent edge.
        self.Pred = None

#---------------------------------------------------------------------------------------------------
# Maintenance.
#---------------------------------------------------------------------------------------------------

    def __enter__(self):
        """
        Function for 'with ... as' context.
        """

        return self

#---------------------------------------------------------------------------------------------------

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Method for 'with ... as' context.

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
            tp -- value of the field 'type'.

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
            nm -- value of the field 'name'.

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
# Parent and children.
#---------------------------------------------------------------------------------------------------

    def IsRoot(self):
        """
        Root check.

        Result:
            True -- if is root,
            False -- if is not root.
        """

        return self.Pred == None

#---------------------------------------------------------------------------------------------------

    def IsLeaf(self):
        """
        List check.

        Result:
            True -- if is leaf,
            False -- if is not leaf.
        """

        return self.Succs == []

#---------------------------------------------------------------------------------------------------

    def ChildrenCount(self):
        """
        Count of children.

        Result:
            Count of children.
        """

        return len(self.Succs)

#---------------------------------------------------------------------------------------------------

    def Parent(self):
        """
        Parent tree.

        Result:
            Parent tree or None.
        """

        if self.IsRoot():
            return None
        else:
            return self.Pred.Parent

#---------------------------------------------------------------------------------------------------

    def Children(self):
        """
        Get all children.

        Result:
            Children.
        """

        return [succ.Child for succ in self.Succs]

#---------------------------------------------------------------------------------------------------

    def Child(self, i):
        """
        Child with given number.

        Arguments:
            i -- number of the child.

        Result:
            Child tree or None.
        """

        if (i >= 0) and (i < self.ChildrenCount()):
            return self.Succs[i].Child
        else:
            return None

#---------------------------------------------------------------------------------------------------

    def FirstChild(self):
        """
        Get first child.

        Result:
            First child tree or None.
        """

        return self.Child(0)

#---------------------------------------------------------------------------------------------------

    def LastChild(self):
        """
        Get last child.

        Result:
            Last child tree or None.
        """

        return self.Child(self.ChildrenCount() - 1)

#---------------------------------------------------------------------------------------------------

    def IsFirstChild(self):
        """
        Check if the three is first child.

        Result:
            True -- if it is the first child,
            False -- otherwise.
        """

        # False for root.
        if self.IsRoot():
            return False;

        # Check.
        return self == self.Parent().FirstChild()

#---------------------------------------------------------------------------------------------------

    def IsLastChild(self):
        """
        Check if the three is lsat child.

        Result:
            True -- if it is the last child,
            False -- otherwise.
        """

        # False for root.
        if self.IsRoot():
            return False;

        # Check.
        return self == self.Parent().LastChild()

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
            return 1 + self.Parent().Level()

#---------------------------------------------------------------------------------------------------
# Elements management.
#---------------------------------------------------------------------------------------------------

    def AddChild(self, t):
        """
        Add new child tree.

        Arguments:
            t -- child tree.

        Result:
            Added child.
        """

        # Links.
        edge = Edge(self, t)

        self.Succs.append(edge);
        t.Pred = edge;

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

        t = Tree(tp, nm, descr);

        return self.AddChild(t);

#---------------------------------------------------------------------------------------------------
# Print.
#---------------------------------------------------------------------------------------------------

    def BaseStr(self):
        """
        Basic string.

        Result:
            Basic string.
        """

        s = '[' + self.Type + ':' + self.Name + ']'

        if self.Descr != '':
            s = s + ' ' + self.Descr

        return s

#---------------------------------------------------------------------------------------------------

    def BaseStrWithEdgePrefix(self):
        """
        Base string with edge prefix.

        Result:
            String.
        """

        if self.IsRoot():
            pref = ''
        else:
            pred = self.Pred
            pref = pred.ToString()

        bs = self.BaseStr()

        if pref == '':
            return bs
        else:
            return pref + ' ' + bs

#---------------------------------------------------------------------------------------------------

    def PropertiesStrings(self):
        """
        Strings of all properties.

        Result:
            Properties string.
        """

        s = []

        # Add all properties.
        for p in self.Dict:
            v = self.Get(p)
            if not (type(v) is type(lambda x: x)):
                s.append(p + ' = ' + str(self.Get(p)))

        return s

#---------------------------------------------------------------------------------------------------

    def FirstIndent(self):
        """
        Get first indentation for print.

        Result:
            First indentation.
        """

        res = ''

        # Init.
        is_arrow = True
        cur = self
        pre = None

        # Main loop.
        while cur != None:

            # Define indentation.
            if is_arrow:
                ind = '--->'
                is_arrow = False
            elif res[:4] == '--->':
                ind = '   |'
            elif (pre != None) and pre.IsLastChild():
                ind = '    '
            else:
                ind = '   |'

            # Grow.
            res = ind + res

            # Shift.
            pre = cur
            cur = cur.Parent()

        res = res[3:] + ' '

        if self.IsRoot():
            res = '#' + res[1:]

        return res

#---------------------------------------------------------------------------------------------------

    def SecondIndent(self):
        """
        Get second indentation for print.

        Result:
            Second indentation.
        """

        res = ""

        # Init.
        cur = self
        pre = None

        # Main loop.
        while cur != None:

            # Define indentation.
            if cur.ChildrenCount() == 0:
                ind = '    '
            elif (pre != None) and pre.IsLastChild():
                ind = '    '
            else:
                ind = '   |'

            # Grow.
            res = ind + res

            # Shift.
            pre = cur
            cur = cur.Parent()

        res = res[3:] + ' '

        return res

#---------------------------------------------------------------------------------------------------

    def Print(self, level, is_recursive):
        """
        Print.

        Arguments:
            level -- tree level,
            is_recursive -- recursive print is needed or not.
        """

        # Properties.
        props_strs = self.PropertiesStrings()
        cnt = len(props_strs)

        # Get indent strings.
        extra_base_str = self.BaseStrWithEdgePrefix()
        if cnt > 0:
            extra_base_str = extra_base_str + ' ('
        space_base_str = ' ' * len(extra_base_str)

        # Indentation.
        f_ind = self.FirstIndent()
        s_ind = self.SecondIndent()

        # Print.
        if cnt == 0:
            print(f_ind + extra_base_str)
        elif cnt == 1:
            print(f_ind + extra_base_str + props_strs[0] + ')')
        else:
            print(f_ind + extra_base_str + props_strs[0])
            for prop_str in props_strs[1:-1]:
                print(s_ind + space_base_str + prop_str)
            print(s_ind + space_base_str + props_strs[cnt - 1] + ')')

        # Print all children.
        if is_recursive:
            for ch in self.Children():
                ch.Print(level + 1, True)

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
            for c in self.Children():
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
# Funcional.
#---------------------------------------------------------------------------------------------------

    def Apply(self, apply_fun, filter_fun = None, is_downward = True):
        """
        Apply 'apply_fun' function to all elements of the tree.

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
            for ch in self.Children():
                ch.Apply(apply_fun, filter_fun, is_downward)
        else:
            # Children.
            for ch in self.Children():
                ch.Apply(apply_fun, filter_fun, is_downward)

            # Apply.
            if is_apply:
                apply_fun(self)

#---------------------------------------------------------------------------------------------------

    def ApplyUpward(self, apply_fun, filter_fun = None):
        """
        Apply 'apply_fun' function to all elements of the tree.

        Arguments:
            apply_fun -- apply function,
            filter_fun -- additional filter function.
        """

        self.Apply(apply_fun, filter_fun, False)

#---------------------------------------------------------------------------------------------------

    def ApplyDownward(self, apply_fun, filter_fun = None):
        """
        Apply 'apply_fun' function to all elements of the tree.

        Arguments:
            apply_fun -- apply function,
            filter_fun -- additional filter function.
        """

        self.Apply(apply_fun, filter_fun, True)

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

            vs = [succ.Child.Get(prop) * succ.GetWithAlternate('count', 1.0)
                  for succ in self.Succs]

            self.Set(prop, sum(vs))

#---------------------------------------------------------------------------------------------------
# Tests.
#---------------------------------------------------------------------------------------------------

if __name__ == '__main__':

    test_number = 1

    if test_number == 1:
        car = Tree('car', 'Toyota', 'Land Cruiser')
        car.Set('max_speed', 300.0)
        car.Set('weight', 2000.0)
        car.Set('color', 'white')
        car.AddChildTN('wheel', 'Bridgestone')
        with car.FindElementByTypeName('wheel', 'Bridgestone') as wheel:
            wheel.Set('color', 'black')
            wheel.Set('height', 5.0)
            wheel.Set('radius', 40.0)
            wheel.InEdge().Set('count', 4)
        car.AddChildTN('engine', 'V8')
        car.PrintTree()
    else:
        pass

#---------------------------------------------------------------------------------------------------
