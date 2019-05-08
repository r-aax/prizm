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

        # X coordinate for visualization.
        self.X = None

        # Node number.
        self.N = None

        # Mark for additional purposes.
        self.Mark = False

        # Data.
        self.Data = data

        # Links to parent and children.
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

    def TreeSizes(self, deltas, margins):
        """
        Tree sizes.

        Arguments:
            deltas -- distances between nodes,
            margins -- margins.

        Result:
            Tree sizes.
        """

        (dx, dy) = deltas
        (mx, my) = margins

        return (int((self.Width() - 1) * dx + 2 * mx), int((self.Height() - 1) * dy + 2 * my))

#---------------------------------------------------------------------------------------------------

    def NodeCoordinates(self, deltas, margins):
        """
        Coordinates of node.

        Arguments:
            deltas -- distances between nodes,
            margins -- margins.

        Result:
            Node coordinates.
        """

        (dx, dy) = deltas
        (mx, my) = margins

        return (int(self.X * dx + mx), int((self.Height() - 1) * dy + my))

#---------------------------------------------------------------------------------------------------
# Clustering.
#---------------------------------------------------------------------------------------------------

def ierarchical_clustering(ps, k = 1):
    """
    Ierarchical clustering.

    Arguments:
        ps -- array of points,
        k -- clusters count for mark.

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

def expand_to_circle(c, r):
    """
    Expand center coordinates to circle corners.

    Arguments:
        c -- coordinates,
        r -- radius.

    Result:
        Expanded coordinates.
    """

    (cx, cy) = c

    return (cx - r, cy - r, cx + r, cy + r)

#---------------------------------------------------------------------------------------------------

def draw_ierarchical_tree(it, deltas = (40, 40), margins = (10, 10)):
    """
    Draw ierarchical tree.

    Arguments:
        it -- ierarchical tree,
        deltas -- distances between nodes,
        margins -- margins.
    """

    # Create image.
    img = Image.new('RGB', it.TreeSizes(deltas, margins), color = (230, 230, 230))
    c = aggdraw.Draw(img)
    c.setantialias(True)

    # Recursive draw.
    draw_ierarchical_tree_on_img(it, c, deltas, margins)

    # Flush, save and show.
    c.flush()
    img.save('itree.png')
    img.show()

#---------------------------------------------------------------------------------------------------

def draw_ierarchical_tree_on_img(it, c, deltas, margins):
    """
    Draw ierarchical tree on image.

    Arguments:
        it -- ierarchical tree,
        c -- canvas,
        deltas -- distances between nodes,
        margins -- margins.
    """

    # Pens and brushes.
    op = aggdraw.Pen('orange', 2.0)
    rp = aggdraw.Pen('red', 2.0)
    bb = aggdraw.Brush('blue')

    # Coordinates.
    p = it.NodeCoordinates(deltas, margins)

    # Draw children.
    for ch in it.Children:
        c.line(p + ch.NodeCoordinates(deltas, margins), op)
        draw_ierarchical_tree_on_img(ch, c, deltas, margins)

    # Draw point.
    c.ellipse(expand_to_circle(p, 3), op, bb)
    if it.Mark:
        c.ellipse(expand_to_circle(p, 10), rp)

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
    draw_ierarchical_tree(tree, deltas = (10, 40), margins = (12, 12))

#---------------------------------------------------------------------------------------------------
