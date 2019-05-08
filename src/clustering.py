# -*- coding: utf-8 -*-
"""
Classes and functions for clustering.

Created on Tue May  7 12:35:16 2019

@author: Rybakov
"""

import numbers
import random
import numpy as np
import aggdraw
from PIL import Image

#---------------------------------------------------------------------------------------------------
# Ierarchical tree class.
#---------------------------------------------------------------------------------------------------

class ITree:
    """
    Ierarchical tree for clustering analysis.
    """

#---------------------------------------------------------------------------------------------------

    def __init__(self, data = None):
        """
        Constructor.

        Arguments:
            data -- data.
        """

        self.X = None
        self.N = None
        self.Mark = False
        self.Data = data
        self.Parent = None
        self.Children= []

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

    def ChildrenCount(self):
        """
        Get children count.

        Result:
            Children count.
        """

        return len(self.Children)

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
            hs = [ch.Height() for ch in self.Children]
            return 1 + max(hs)

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
            ws = [ch.Width() for ch in self.Children]
            return sum(ws)

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
            RIght leaf.
        """

        if self.IsLeaf():
            return self
        else:
            return self.RightChild().RightLeaf()

#---------------------------------------------------------------------------------------------------

    def Str(self):
        """
        String representation.

        Result:
            String.
        """

        if self.Data == None:
            return 'None'
        elif isinstance(self.Data, numbers.Number):
            return 'x = %d, n = %d, data = %f (1d)' % (self.X, self.N, self.Data)
        else:
            return 'str'

#---------------------------------------------------------------------------------------------------

    def Print(self):
        """
        Print tree.
        """

        level = self.Level()
        indent = ' ' * level
        print('%s [L%d] : %s' % (indent, level, self.Str()))

        # Print all children.
        for c in self.Children:
            c.Print()

#---------------------------------------------------------------------------------------------------
# Clustering.
#---------------------------------------------------------------------------------------------------

def ierarchical_clustering(ps, k = 1):
    """
    Ierarchical clustering.

    Arguments:
        ps -- array of points,
        k -- clusters count.

    Result:
        Clustering tree.
    """

    trees = [ITree(p) for p in ps]

    # Now set pseudocoordinates and numbers.
    for i in range(len(trees)):
        trees[i].X = i
        trees[i].N = i

    # Next number.
    next_n = len(trees)

    while len(trees) > 1:

        # Find two nearest points.
        fpi = 0
        dist = abs(trees[fpi].Data - trees[fpi + 1].Data)
        for i in range(fpi + 1, len(trees) - 1):
            new_dist = abs(trees[i].Data - trees[i + 1].Data)
            if new_dist < dist:
                fpi = i
                dist = new_dist

        # Merge points first_point_index and first_point_index + 1.
        ch1 = trees[fpi]
        ch2 = trees[fpi + 1]
        chs = [ch1, ch2]
        new_tree = ITree()
        new_tree.Children = chs
        new_tree.Data = 0.5 * (ch1.Data + ch2.Data)
        new_tree.X = 0.5 * (ch1.X + ch2.X)
        new_tree.N = next_n
        next_n = next_n + 1
        ch1.Parent = new_tree
        ch2.Parent = new_tree
        trees = trees[ : fpi] + [new_tree] + trees[fpi + 2 : ]

        # Mark.
        if len(trees) == k:
            for tree in trees:
                tree.Mark = True

    # Last tree is a root.
    return trees[0]

#---------------------------------------------------------------------------------------------------
# Visualization.
#---------------------------------------------------------------------------------------------------

def draw_ierarchical_tree(it, dx = 40, dy = 40, mx = 10, my = 10):
    """
    Draw ierarchical tree.

    Arguments:
        it -- ierarchical tree,
        dx -- horizontal distance between two neighbours,
        dy -- vertical distance between nodes,
        mx -- horizontal margin,
        my -- vertical margin.
    """

    # Create image.
    width = int((it.Width() - 1) * dx + 2 * mx)
    height = int((it.Height() - 1) * dy + 2 * my)
    img = Image.new('RGB', (width, height), color = (230, 230, 230))
    c = aggdraw.Draw(img)
    c.setantialias(True)

    # Recursive draw.
    draw_ierarchical_tree_on_img(it, c, dx, dy, mx, my)

    # Flush, save and show.
    c.flush()
    img.save('itree.png')
    img.show()

#---------------------------------------------------------------------------------------------------

def draw_ierarchical_tree_on_img(it, c, dx, dy, mx, my):
    """
    Draw ierarchical tree on image.

    Arguments:
        it -- ierarchical tree,
        c -- canvas,
        dx -- horizontal distance between two neighbours,
        dy -- vertical distance between nodes,
        mx -- horizontal margin,
        my -- vertical margin.
    """

    # Coordinates.
    cx = int(it.X * dx + mx)
    cy = int((it.Height() - 1) * dy + my)

    # Draw children.
    for ch in it.Children:
        chx = int(ch.X * dx + mx)
        chy = int((ch.Height() - 1) * dy + my)
        c.line((cx, cy, chx, chy), aggdraw.Pen('orange', 2.0))
        draw_ierarchical_tree_on_img(ch, c, dx, dy, mx, my)

    # Draw point.
    c.ellipse((cx - 3, cy - 3, cx + 3, cy + 3),
              aggdraw.Pen('orange', 2.0), aggdraw.Brush('blue'))
    if it.Mark:
        c.ellipse((cx - 10, cy - 10, cx + 10, cy + 10),
                  aggdraw.Pen('red', 2.0))

#---------------------------------------------------------------------------------------------------
# Tests.
#---------------------------------------------------------------------------------------------------

if __name__ == '__main__':

    # Simple range.
    #ps = list(range(150))

    # Random numbers.
    #ps = [100.0 * random.random() for j in range(100)]

    # Normal distribution.
    #ps = [random.normalvariate(50.0, 20.0) for j in range(100)]

    # Normal distribution.
    ps = [random.gammavariate(50.0, 20.0) for j in range(150)]

    ps.sort()
    tree = ierarchical_clustering(ps, k = 12)
    tree.Print()
    draw_ierarchical_tree(tree, dx = 10, mx = 12, my = 12)

#---------------------------------------------------------------------------------------------------
