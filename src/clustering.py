# -*- coding: utf-8 -*-
"""
Classes and functions for clustering.

Created on Tue May  7 12:35:16 2019

@author: Rybakov
"""

import numbers
import math
import random
import numpy as np
import aggdraw
import itertools
from enum import Enum
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

        # Mark for additional purposes.
        self.Mark = False

        # Cluster number.
        self.KN = -1

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
            return 'x = %d, kn = %d, data = %f (1d)' % (self.X, self.KN, self.Data)
        elif isinstance(self.Data, tuple):
            if len(self.Data) == 2:
                return 'x = %d, kn = %d, data = (%f, %f) (2d)' \
                       % (self.X, self.KN, self.Data[0], self.Data[1])
            else:
                raise Exception('wrong data for string representation')
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

        dx, dy = deltas
        mx, my = margins

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

        dx, dy = deltas
        mx, my = margins

        return (int(self.X * dx + mx), int((self.Height() - 1) * dy + my))

#---------------------------------------------------------------------------------------------------

    def Merge(t1, t2):
        """
        Merge two trees.

        Arguments:
            t1 -- first tree,
            t2 -- second tree.

        Result:
            New tree.
        """

        # New tree.
        t = ITree()

        # Links.
        t.Children = [t1, t2]
        t1.Parent = t
        t2.Parent = t

        # Weights and mean function.
        w1 = t1.Width()
        w2 = t2.Width()
        m = lambda v1, v2: (v1 * w1 + v2 * w2) / (w1 + w2)

        # Data.
        if isinstance(t1.Data, numbers.Number):
            t.Data = m(t1.Data, t2.Data)
        elif isinstance(t1.Data, tuple):
            if len(t1.Data) == 2:
                t.Data = (m(t1.Data[0], t2.Data[0]), m(t1.Data[1], t2.Data[1]))
            else:
                raise Exception('wrong tuple len for data while merging')
        else:
            raise Exception('wrong data type for merge')

        # Coordinate X scaled by Width.
        t.X = m(t1.X, t2.X)

        return t

#---------------------------------------------------------------------------------------------------
# Clustering visualization type.
#---------------------------------------------------------------------------------------------------

class ClusteringDrawingType(Enum):

    # Connect nodes with simple lines.
    Lines = 1

    # Connect nodes with orthogonal lines.
    Orthogonal = 2

#---------------------------------------------------------------------------------------------------
# Clustering.
#---------------------------------------------------------------------------------------------------

def metric_data_dist(a, b, r):
    """
    Metric - distance between points.

    Arguments:
        a -- first data,
        b -- second data,
        r -- power.

    Result:
        Metric result.
    """

    if isinstance(a, numbers.Number):
        return abs(a - b)
    elif isinstance(a, tuple):
        return pow(sum([pow(abs(ai - bi), r) for (ai, bi) in zip(a, b)]), 1.0 / r)
    else:
        raise Exception('wrong data for metric_data_dist')

#---------------------------------------------------------------------------------------------------

def metric_tree_dist(t1, t2, r):
    """
    Metric - distance between trees.

    Arguments:
        t1 -- first tree,
        t2 -- second tree,
        r -- power.

    Result:
        Metric result.
    """

    return metric_data_dist(t1.Data, t2.Data, r)

#---------------------------------------------------------------------------------------------------

def ierarchical_clustering(ps,
                           k = 1,
                           metric = lambda t1, t2: metric_tree_dist(t1, t2, 1.0)):
    """
    Ierarchical clustering.

    Arguments:
        ps -- array of points,
        k -- clusters count for mark,
        metric -- metric function.

    Result:
        Clustering tree.
    """

    trees = [ITree(p) for p in ps]

    # Now set pseudocoordinates and numbers.
    for i in range(len(trees)):
        trees[i].X = i
        trees[i].NN = i

    # Next number.
    next_n = len(trees)

    while len(trees) > 1:

        # Find two nearest points.
        fpi = 0
        m = metric(trees[fpi], trees[fpi + 1])
        for i in range(1, len(trees) - 1):
            new_m = metric(trees[i], trees[i + 1])
            if new_m < m:
                fpi, m = i, new_m

        # Merge points first_point_index and first_point_index + 1.
        ch1, ch2 = trees[fpi], trees[fpi + 1]
        new_tree = ITree.Merge(ch1, ch2)
        new_tree.N = next_n
        next_n = next_n + 1
        trees = trees[ : fpi] + [new_tree] + trees[fpi + 2 : ]

        # Mark.
        if len(trees) == k:
            for i in range(len(trees)):
                trees[i].KN = i
                trees[i].Mark = True

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

def random_color():
    """
    Random 3-component color.

    Result:
        Random color.
    """

    return (random.randint(0, 255), random.randint(0, 255), random.randint(0, 255))

#---------------------------------------------------------------------------------------------------

def pretty_color(n):
    """
    Pretty color by number.

    Arguments:
        n -- number.

    Result:
        Color.
    """

    colors = [(146,  18,  80), ( 21, 161,  44), ( 51,  83, 142), (178, 143, 131),
              (227,  50,  11), (113, 144, 184), ( 12,  75,  20), ( 71, 103, 220),
              (220, 129, 118), (126, 110,  20), (199, 114, 146), ( 34,  35, 177)]
    colors_count = len(colors)

    if n < colors_count:
        pc = colors[n]
    else:
        pc = random_color()

    return pc

#---------------------------------------------------------------------------------------------------

def draw_trajectory(ps):
    """
    Draw trajectory.

    Arguments:
        ps -- points array.
    """

    # Transform image.
    w = 600
    h = 200
    sw = 100
    sh = 80

    # Create image.
    img = Image.new('RGB', (w, h), color = (230, 230, 230))
    c = aggdraw.Draw(img)
    c.setantialias(True)

    # Draw points.
    for i in range(len(ps) - 1):
        p1 = ps[i]
        p2 = ps[i + 1]
        c.line((p1[0] * sw,
                p1[1] * sh + h / 2,
                p2[0] * sw,
                p2[1] * sh + h / 2),
               aggdraw.Pen('red', 1.0))

    # Flush save and show.
    c.flush()
    img.save('trajectory.png')
    img.show()

#---------------------------------------------------------------------------------------------------


def draw_ierarchical_tree(it,
                          deltas = (40, 40),
                          margins = (10, 10),
                          pen = aggdraw.Pen('orange', 2.0),
                          drawing_type = ClusteringDrawingType.Lines):
    """
    Draw ierarchical tree.

    Arguments:
        it -- ierarchical tree,
        deltas -- distances between nodes,
        margins -- margins,
        pen -- pen,
        drawing_type -- drawing type.
    """

    # Create image.
    img = Image.new('RGB', it.TreeSizes(deltas, margins), color = (230, 230, 230))
    c = aggdraw.Draw(img)
    c.setantialias(True)

    # Recursive draw.
    draw_ierarchical_tree_on_img(it, c, deltas, margins, pen, drawing_type)

    # Flush, save and show.
    c.flush()
    img.save('itree.png')
    img.show()

#---------------------------------------------------------------------------------------------------

def draw_ierarchical_tree_on_img(it, c, deltas, margins, pen, drawing_type):
    """
    Draw ierarchical tree on image.

    Arguments:
        it -- ierarchical tree,
        c -- canvas,
        deltas -- distances between nodes,
        margins -- margins,
        pen -- pen,
        drawing_type -- drawing type.
    """

    # Pens and brushes.
    mark_pen = aggdraw.Pen('red', 2.0)
    brush = aggdraw.Brush('silver')

    # Change color for subtree.
    if it.Mark:
        pen = aggdraw.Pen(pretty_color(it.KN), 2.0)

    # Coordinates.
    p = it.NodeCoordinates(deltas, margins)

    # Draw children.
    for ch in it.Children:
        chp = ch.NodeCoordinates(deltas, margins)

        if drawing_type == ClusteringDrawingType.Lines:
            c.line(p + chp, pen)
        elif drawing_type == ClusteringDrawingType.Orthogonal:
            (px, py) = p
            (chpx, chpy) = chp
            c.line((px, py, chpx, py, chpx, chpy), pen)
        else:
            raise Exception('wrong drawing type : %s' % str(drawing_type))

        draw_ierarchical_tree_on_img(ch, c, deltas, margins, pen, drawing_type)

    # Draw point.
    c.ellipse(expand_to_circle(p, 3), pen, brush)
    if it.Mark:
        c.ellipse(expand_to_circle(p, 10), mark_pen)

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

    # Gamma distribution.
    #ps = [random.gammavariate(50.0, 20.0) for j in range(150)]

    ps = [(x / 100.0, math.sin(x / 25.0)) for x in range(600)]
    ps[300] = (ps[300][0], 1000.0)
    ps.sort()
    clust = True

    # Just draw or clluster.
    if not clust:
        draw_trajectory(ps)
    else:
        tree = ierarchical_clustering(ps, k = 12)
        tree.Print()
        draw_ierarchical_tree(tree, deltas = (10, 40), margins = (12, 12),
                              drawing_type = ClusteringDrawingType.Orthogonal)

#---------------------------------------------------------------------------------------------------
