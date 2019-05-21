# -*- coding: utf-8 -*-
"""
Hierarchical tree.

Created on Tue May 21 16:40:37 2019

@author: Rybakov
"""

import operator
from functools import reduce

class HTree:
    """
    Hierarchical tree.
    """

#---------------------------------------------------------------------------------------------------
# Constructor.
#---------------------------------------------------------------------------------------------------

    def __init__(self):
        """
        Constructor.
        """

        # Links to parent and children.
        self.Parent = None
        self.Children= []

#---------------------------------------------------------------------------------------------------
# Basic functions.
#---------------------------------------------------------------------------------------------------

    def ChildrenCount(self):
        """
        Get children count.

        Result:
            Children count.
        """

        return len(self.Children)

#---------------------------------------------------------------------------------------------------

    def IsRoot(self):
        """
        Root check.

        Result:
            True -- if is a root,
            False -- if is not a root.
        """

        return self.Parent == None

#---------------------------------------------------------------------------------------------------

    def IsLeaf(self):
        """
        Leaf check.

        Result:
            True -- if is a leaf,
            False -- if is not a leaf.
        """

        return self.Children == []

#---------------------------------------------------------------------------------------------------

    def Level(self):
        """
        Level.

        Result:
            Level.
        """

        if self.IsRoot():
            return 0
        else:
            return 1 + self.Parent.Level()

#---------------------------------------------------------------------------------------------------

    def Height(self):
        """
        Tree height (levels count).

        Result:
            Height.
        """

        if self.IsLeaf():
            return 1
        else:
            return 1 + max([ch.Height() for ch in self.Children])

#---------------------------------------------------------------------------------------------------

    def Width(self):
        """
        Tree width (leafs count).

        Result:
            Width.
        """

        if self.IsLeaf():
            return 1
        else:
            return sum([ch.Width() for ch in self.Children])

#---------------------------------------------------------------------------------------------------

    def LeftChild(self):
        """
        Get left child.

        Result:
            Left child.
        """

        if self.IsLeaf():
            return None
        else:
            return self.Children[0]

#---------------------------------------------------------------------------------------------------

    def RightChild(self):
        """
        Get right child.

        Result:
            Right child.
        """

        if self.IsLeaf():
            return None
        else:
            return self.Children[self.ChildrenCount() - 1]

#---------------------------------------------------------------------------------------------------

    def LeftLeaf(self):
        """
        Get left leaf.

        Result:
            Left leaf.
        """

        if self.IsLeaf():
            return self
        else:
            return self.LeftChild().LeftLeaf()

#---------------------------------------------------------------------------------------------------

    def RightLeaf(self):
        """
        Get right leaf.

        Result:
            Right leaf.
        """

        if self.IsLeaf():
            return self
        else:
            return self.RightChild.RightLeaf()

#---------------------------------------------------------------------------------------------------

    def NextLeafLeftRoRight(self, leaf):
        """
        Get next leaf (left to right walk).

        Arguments:
            leaf -- leaf.

        Result:
            Next leaf or none.
        """

        if not leaf.IsLeaf():
            raise Exception('not a leaf')

        # Find next leaf.
        cur = leaf
        while True:

            if cur.IsRoot():
                return None

            p = cur.Parent

            if cur == p.RightChild():
                cur = p
            else:
                for i in range(p.ChildrenCount()):
                    if cur == p.Children[i]:
                        return p.Children[i + 1].LeftLeaf()

#---------------------------------------------------------------------------------------------------
# Functional.
#---------------------------------------------------------------------------------------------------

    def Apply(self, fun):
        """
        Apply function to tree recursive.

        Arguments:
            fun -- function.
        """

        fun(self)
        for ch in self.Children:
            ch.Apply(fun)

#---------------------------------------------------------------------------------------------------

    def MapLeafs(self, fun):
        """
        Apply function to all leafs and return result in list.

        Arguments:
            fun -- map function.

        Result:
            Map result.
        """

        if self.IsLeaf():
            return fun(self)
        else:
            return reduce(operator.__concat__,
                          [ch.MapLeafs(fun) for ch in self.Children],
                          [])

#---------------------------------------------------------------------------------------------------
