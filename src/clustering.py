# -*- coding: utf-8 -*-
"""
Classes and functions for clustering.

Created on Tue May  7 12:35:16 2019

@author: Rybakov
"""

import numbers
import math
import fun
import random
import numpy as np
import aggdraw
import itertools
import mth
from functools import reduce
from draw import Drawer
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

        # Height difference.
        self.HeightDiff = None

        # Is overshoot.
        self.IsOvershoot = False

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

        height_diff_str = (', hd = ' + str(self.HeightDiff)) if self.HeightDiff != None else ''

        if self.Data == None:
            return 'None'
        elif isinstance(self.Data, numbers.Number):
            return 'x = %f, kn = %d, data = %f%s (1d)' \
                   % (self.X, self.KN, self.Data, height_diff_str)
        elif isinstance(self.Data, tuple):
            if len(self.Data) == 2:
                return 'x = %f, kn = %d, data = (%f, %f)%s (2d)' \
                       % (self.X, self.KN, self.Data[0], self.Data[1], height_diff_str)
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
        ws = [t1.Width(), t2.Width()]

        # Data.
        if isinstance(t1.Data, numbers.Number):
            t.Data = mth.avg_weighted([t1.Data, t2.Data], ws)
        elif isinstance(t1.Data, tuple):
            if len(t1.Data) == 2:
                t.Data = (mth.avg_weighted([t1.Data[0], t2.Data[0]], ws),
                          mth.avg_weighted([t1.Data[1], t2.Data[1]], ws))
            else:
                raise Exception('wrong tuple len for data while merging')
        else:
            raise Exception('wrong data type for merge')

        # Coordinate X scaled by Width.
        t.X = mth.avg_weighted([t1.X, t2.X], ws)

        return t

#---------------------------------------------------------------------------------------------------

    def CalculateHeightDifferences(self):
        """
        Calculate height differences.
        """

        if self.IsLeaf():
            self.HeightDiff = self.Parent.Height() - self.Height()
        else:
            for ch in self.Children:
                ch.CalculateHeightDifferences()

#---------------------------------------------------------------------------------------------------

    def MaxOvershootLeaf(self):
        """
        Maximum overshoot leaf.

        Result:
            Maximum overshoot leaf.
        """

        if self.IsLeaf():
            if self.IsOvershoot:
                return None
            else:
                return self
        else:
            ms = [ch.MaxOvershootLeaf() for ch in self.Children]
            cur_m = None
            for m in ms:
                if cur_m == None:
                    cur_m = m
                elif m == None:
                    pass
                elif m.HeightDiff > cur_m.HeightDiff:
                    cur_m = m
                else:
                    pass

        return cur_m

#---------------------------------------------------------------------------------------------------

    def FindOvershoots(self, n):
        """
        Find overshoots.

        Arguments:
            n -- count.
        """

        for i in range(n):
            m = self.MaxOvershootLeaf()
            m.IsOvershoot = True

#---------------------------------------------------------------------------------------------------

    def SetLeafsXs(self, start = 0):
        """
        Set X coordinates for leafs.

        Arguments:
            start -- start position
        """

        if self.IsLeaf():
            self.X = start
        else:
            cur = start
            for i in range(self.ChildrenCount()):
                ch = self.Children[i]
                ch.SetLeafsXs(cur)
                cur = cur + ch.Width()

#---------------------------------------------------------------------------------------------------

    def RefreshXs(self):
        """
        Refresh Xs.
        """

        if self.IsLeaf():
            pass
        else:
            for ch in self.Children:
                ch.RefreshXs()
            xws = [(ch.X, ch.Width()) for ch in self.Children]
            (xs, ws) = fun.unzip(xws)
            self.X = mth.avg_weighted(xs, ws)

#---------------------------------------------------------------------------------------------------
# Clustering nearest finding type.
#---------------------------------------------------------------------------------------------------

class ClusteringNearestType(Enum):

    # Find only for adjacent elements (good for 1d clustering).
    OnlyAdjacent = 1

    # All elements.
    All = 2

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

def metric_data_idist(a, b, i):
    """
    Metric - distance between i-th components.

    Arguments:
        a -- first data,
        b -- second data,
        i -- element number.

    Result:
        Metric result.
    """

    if not isinstance(a, tuple):
        raise Exception('wrong data for metric_data_idist')

    return abs(a[i] - b[i])

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

def metric_tree_idist(t1, t2, i):
    """
    Metric - distance between i-th components.

    Arguments:
        t1 -- first tree,
        t2 -- second tree,
        i -- element number.

    Result:
        Metric result.
    """

    return metric_data_idist(t1.Data, t2.Data, i)

#---------------------------------------------------------------------------------------------------

def ierarchical_clustering(ps,
                           k = 1,
                           metric = lambda t1, t2: metric_tree_dist(t1, t2, 1.0),
                           nearest_type = ClusteringNearestType.All):
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

        # Initial indexes and metric.
        fpi, spi, m = 0, 1, metric(trees[0], trees[1])

        # Find two nearest points.
        if nearest_type == ClusteringNearestType.OnlyAdjacent:
            for i in range(1, len(trees) - 1):
                new_m = metric(trees[i], trees[i + 1])
                if new_m < m:
                    fpi, spi, m = i, i + 1, new_m
        elif nearest_type == ClusteringNearestType.All:
            for i in range(0, len(trees)):
                for j in range(i + 1, len(trees)):
                    new_m = metric(trees[i], trees[j])
                    if new_m < m:
                        fpi, spi, m = i, j, new_m
        else:
            raise Exception('wrong nearest type for clustering')

        # Merge two trees.
        new_tree = ITree.Merge(trees[fpi], trees[spi])
        new_tree.N = next_n
        next_n = next_n + 1
        trees = trees[ : fpi] + [new_tree] + trees[fpi + 1 : spi] + trees[spi + 1 : ]

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

def draw_data(ps,
              draw_trajectory = False,
              clusters = None,
              pic_size = (640, 480),
              is_axis = True,
              grid = None,
              filename = None):
    """
    Draw data.

    Arguments:
        ps -- points array,
        pic_sizr -- picture size,
        is_axis -- need to draw axis.
    """

    # Points characteristics.
    min_coords = reduce(lambda p1, p2: (min(p1[0], p2[0]), min(p1[1], p2[1])), ps[1:], ps[0])
    max_coords = reduce(lambda p1, p2: (max(p1[0], p2[0]), max(p1[1], p2[1])), ps[1:], ps[0])

    # Drawer ini.
    D = Drawer(draw_area = min_coords + max_coords, pic_size = pic_size)

    # Axis.
    if is_axis:
        D.Axis()

    # Grid.
    if grid != None:
        (min_x, min_y) = min_coords
        (max_x, max_y) = max_coords
        pen = aggdraw.Pen('silver', 1.0)
        (gx, gy) = grid
        gx_cur = gx
        while gx_cur <= max_x:
            D.Line((gx_cur, min_y), (gx_cur, max_y), pen = pen)
            gx_cur = gx_cur + gx
        gx_cur = -gx
        while gx_cur >= min_x:
            D.Line((gx_cur, min_y), (gx_cur, max_y), pen = pen)
            gx_cur = gx_cur - gx
        gy_cur = gy
        while gy_cur <= max_y:
            D.Line((min_x, gy_cur), (max_x, gy_cur), pen = pen)
            gy_cur = gy_cur + gy
        gy_cur = -gy
        while gy_cur >= min_y:
            D.Line((min_x, gy_cur), (max_x, gy_cur), pen = pen)
            gy_cur = gy_cur - gy

    # Draw points.
    pen = aggdraw.Pen('red', 1.0)
    point_pen = pen
    point_brush = aggdraw.Brush('red')
    if draw_trajectory:
        for i in range(len(ps) - 1):
            D.Line(ps[i], ps[i + 1], pen = pen)
    for i in range(len(ps)):
        if clusters != None:
            color = pretty_color(clusters[i])
            point_pen = aggdraw.Pen(color, 1.0)
            point_brush = aggdraw.Brush(color)
        D.Point(ps[i], 3, pen = point_pen, brush = point_brush)

    # Flush save and show.
    D.FSS(filename = filename)

#---------------------------------------------------------------------------------------------------


def draw_ierarchical_tree(it,
                          deltas = (40, 40),
                          margins = (10, 10),
                          pen = aggdraw.Pen('orange', 2.0),
                          drawing_type = ClusteringDrawingType.Lines,
                          filename = None):
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
    if filename != None:
        img.save(filename)
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
    if it.IsOvershoot:
        c.ellipse(expand_to_circle(p, 5), aggdraw.Pen('black', 2.0), aggdraw.Brush('black'))
    else:
        c.ellipse(expand_to_circle(p, 3), pen, brush)
    if it.Mark:
        c.ellipse(expand_to_circle(p, 10), mark_pen)

#---------------------------------------------------------------------------------------------------
# Tests.
#---------------------------------------------------------------------------------------------------

def test_set_points2d(n, k):
    """
    Test set with 2D points.

    Arguments:
        n -- points count,
        k -- clusters count.

    Result:
        Points.
    """

    # Random point.
    rp = lambda : (random.uniform(0.0, 100.0), random.uniform(0.0, 100.0))

    # Clusters.
    ks = [(rp(), random.uniform(0.0, 20.0)) for x in range(k)]

    # Generare points.
    ps = []
    for i in range(n):
        ki = random.randint(0, k - 1)
        ((cx, cy), r) = ks[ki]
        while True:
            (px, py) = (random.uniform(cx - r, cx + r), random.uniform(cy - r, cy + r))
            dx = px - cx
            dy = py - cy
            if dx * dx + dy * dy < r * r:
                break
        ps = ps + [(px, py)]

    return ps

#---------------------------------------------------------------------------------------------------

def test_set_trajectory(ripple, overshoot):
    """
    Test set Trajectory.

    Arguments:
        ripple -- ripple intensive,
        overshoot -- overshoot value.

    Result:
        Test set.
    """

    s = 200
    ps = [(x / 15.0, math.sin(x / 15.0)) for x in range(s)]

    # Ripple.
    for i in range(s):
        dx = random.uniform(-ripple, ripple)
        dy = random.uniform(-ripple, ripple)
        ps[i] = (ps[i][0] + dx, ps[i][1] + dy)

    # Overshoots.
    for i in range(4):
        ind = 40 * (i + 1)
        dx = random.uniform(-overshoot, overshoot)
        dy = random.uniform(-overshoot, overshoot)
        ps[ind] = (ps[ind][0] + dx, ps[ind][1] + dy)

    return ps

#---------------------------------------------------------------------------------------------------
# Run type.
#---------------------------------------------------------------------------------------------------

class RunType(Enum):
    """
    Type of test run.
    """

    # Points 2D.
    Points2D = 1

    # Trajectory test.
    Trajectory = 2

#---------------------------------------------------------------------------------------------------


if __name__ == '__main__':

    points_count = 200
    run = RunType.Points2D

    if run == RunType.Points2D:
        ps = test_set_points2d(points_count, k = 12)
        draw_data(ps,
                  pic_size = (600, 600),
                  grid = (10.0, 10.0),
                  filename = 'points2d_init.png')
        #tree = ierarchical_clustering(ps, k = 12)
        #tree.SetLeafsXs()
        #tree.RefreshXs()
        #tree.CalculateHeightDifferences()
        #tree.FindOvershoots(5)
        #tree.Print()
        #draw_ierarchical_tree(tree, deltas = (10, 40), margins = (12, 12),
        #                      drawing_type = ClusteringDrawingType.Orthogonal,
        #                      filename = 'points2d_tree.png')
    elif run == RunType.Trajectory:
        ps = test_set_trajectory(0.03, 0.9)
        draw_data(ps,
                  draw_trajectory = True,
                  pic_size = (1000, 600),
                  grid = (1.0, 0.3),
                  filename = 'trajectory_init.png')
        tree = ierarchical_clustering(ps, k = 12)
        tree.CalculateHeightDifferences()
        tree.FindOvershoots(4)
        tree.Print()
        draw_ierarchical_tree(tree, deltas = (10, 40), margins = (12, 12),
                              drawing_type = ClusteringDrawingType.Orthogonal,
                              filename = 'trajectory_tree.png')
    else:
        raise Exception('wrong test run type')

#---------------------------------------------------------------------------------------------------
