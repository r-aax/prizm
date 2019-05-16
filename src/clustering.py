# -*- coding: utf-8 -*-
"""
Classes and functions for clustering.

Created on Tue May  7 12:35:16 2019

@author: Rybakov
"""

import numbers
import math
import fun
import lst
import operator
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

    def NextLeafLeftRoRight(self, leaf):
        """
        Get next leaf (left to right walk).

        Arguments:
            leaf - leaf.

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

    def ClusterNumber(self):
        """
        Node cluster number.

        Result:
            Cluster number.
        """

        kn = self.KN

        if self.IsRoot():
            return kn
        elif kn != -1:
            return kn
        else:
            return self.Parent.ClusterNumber()

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

        # Calculate own heights difference.
        if self.IsRoot():
            self.HeightDiff = 0
        elif self.ClusterNumber() == -1:
            self.HeightDiff = 0
        else:
            self.HeightDiff = self.Parent.Height() - self.Height()

        # Calculate childrens' heights differences.
        if not self.IsLeaf():
            for ch in self.Children:
                ch.CalculateHeightDifferences()

#---------------------------------------------------------------------------------------------------

    def MaxOvershootLeaf(self):
        """
        Maximum overshoot leaf.

        Result:
            Maximum overshoot leaf.
        """

        if self.IsOvershoot:
            return None
        elif self.IsLeaf():
            return self
        else:
            ms = [ch.MaxOvershootLeaf() for ch in self.Children]
            cur_m = self
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

    def SetOvershootRecursive(self):
        """"
        Set overshoot recursive.
        """

        self.IsOvershoot = True

        for ch in self.Children:
            ch.SetOvershootRecursive()

#---------------------------------------------------------------------------------------------------

    def FindOvershoots(self, n):
        """
        Find overshoots.

        Arguments:
            n -- count.
        """

        for i in range(n):
            m = self.MaxOvershootLeaf()
            m.SetOvershootRecursive()

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

    def AllData(self):
        """
        All data of tree (leafs).

        Result:
            Data as a list.
        """

        if self.IsLeaf():
            return [self.Data]
        else:
            return reduce(operator.__concat__,
                          [ch.AllData() for ch in self.Children],
                          [])

#---------------------------------------------------------------------------------------------------

    def OvershootData(self):
        """
        Overshootdata (leafs).

        Result:
            Overshoot data.
        """

        if self.IsLeaf():
            if self.IsOvershoot:
                return [self.Data]
            else:
                return []
        else:
            return reduce(operator.__concat__,
                          [ch.OvershootData() for ch in self.Children],
                          [])

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

def metric_data_dist(a, b, r = 2.0):
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

def metric_tree_dist(t1, t2, r = 2.0):
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

def dists_list(t1, t2, r = 2.0):
    """
    List of distances for all pair of leafs.

    Arguments:
        t1 -- first tree,
        t2 -- second tree,
        r -- power.

    Result:
        List of distances between nodes from two trees.
    """

    d1 = t1.AllData()
    d2 = t2.AllData()
    datas = lst.descartes_product(d1, d2)

    return [metric_data_dist(p1, p2, r) for (p1, p2) in datas]

#---------------------------------------------------------------------------------------------------

def metric_tree_min_dist(t1, t2, r = 2.0):
    """
    Metric - minimal distance between nodes.

    Arguments:
        t1 -- first tree,
        t2 -- second tree,
        r -- power.

    Result:
        Metric result.
    """

    return min(dists_list(t1, t2, r))

#---------------------------------------------------------------------------------------------------

def metric_tree_max_dist(t1, t2, r = 2.0):
    """
    Metric - maximal distance between nodes.

    Arguments:
        t1 -- first tree,
        t2 -- second tree,
        r -- power.

    Result:
        Metric result.
    """

    return max(dists_list(t1, t2, r))

#---------------------------------------------------------------------------------------------------

def metric_tree_avg_dist(t1, t2, r = 2.0):
    """
    Metric - average distance between nodes.

    Arguments:
        t1 -- first tree,
        t2 -- second tree,
        r -- power.

    Result:
        Metric result.
    """

    return mth.avg_arith(dists_list(t1, t2, r))

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
                           metric = lambda t1, t2: metric_tree_dist(t1, t2, 2.0),
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

def draw_data(tree,
              draw_trajectory = False,
              draw_clusters = True,
              pic_size = (640, 480),
              is_axis = True,
              grid = None,
              filename = None):
    """
    Draw data.

    Arguments:
        tree -- clustering tree,
        pic_size -- picture size,
        is_axis -- need to draw axi,
        grid -- grid lines characteristics,
        filename -- file name for picture.
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
        D.Grid(grid)

    # Draw clusters lines.
    if draw_clusters:
        leaf1 = tree.LeftLeaf()
        while leaf1 != None:
            leaf2 = tree.NextLeafLeftRoRight(leaf1)
            while leaf2 != None:

                # Draw.
                d1 = leaf1.Data
                d2 = leaf2.Data
                kn1 = leaf1.ClusterNumber()
                kn2 = leaf2.ClusterNumber()
                is_not_overshoots = (not leaf1.IsOvershoot) and (not leaf2.IsOvershoot)
                if is_not_overshoots and (kn1 != -1) and (kn1 == kn2):
                    D.Line(d1, d2, pen = aggdraw.Pen(pretty_color(kn1), 1.0))

                leaf2 = tree.NextLeafLeftRoRight(leaf2)
            leaf1 = tree.NextLeafLeftRoRight(leaf1)

    foreground_color = (230, 230, 230)
    foreground_pen = aggdraw.Pen(foreground_color, 1.0)
    foreground_brush = aggdraw.Brush(foreground_color)

    # Draw points.
    leaf = tree.LeftLeaf()
    while leaf != None:

        # Default colors and etc.
        color = 'red'
        point_radius = 3

        # Define color if cluster number is set.
        if draw_clusters:
            kn = leaf.ClusterNumber()
            if leaf.IsOvershoot:
                color = 'black'
                point_radius = 5
            elif kn != -1:
                color = pretty_color(kn)

        # Define pen, brush and draw point.
        point_pen = aggdraw.Pen(color, 1.0)
        point_brush = aggdraw.Brush(color)
        D.Point(leaf.Data, point_radius, pen = point_pen, brush = point_brush)
        if draw_clusters:
            if not leaf.IsOvershoot:
                D.Point(leaf.Data, 1, pen = foreground_pen, brush = foreground_brush)
        leaf = tree.NextLeafLeftRoRight(leaf)

    # Flush save and show.
    D.FSS(filename = filename)

#---------------------------------------------------------------------------------------------------

def draw_overshoots(ps, overshoots, ok,
                    pic_size = (640, 480),
                    is_axis = True,
                    grid = None,
                    filename = None):
    """
    Draw overshoots.

    Arguments:
        ps -- points,
        overshoots -- list of overshoots,
        ok -- overshoots count (leaders),
        pic_size -- picture size,
        is_axis -- need to draw axis,
        grid -- grid lines characteristics,
        filename -- file name for picture.
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
        D.Grid(grid)

    # Draw points.
    red_pen = aggdraw.Pen('red', 1.0)
    red_brush = aggdraw.Brush('red')
    for p in ps:
        D.Point(p, 3, red_pen, red_brush)

    # Draw overshoots.
    black_pen = aggdraw.Pen('black', 2.0)
    steelblue_brush = aggdraw.Brush('steelblue')
    for overshoot in overshoots[ok:]:
        (r, p) = overshoot
        D.Point(p, 3 * r, black_pen)
    for overshoot in overshoots[:ok]:
        (r, p) = overshoot
        D.Point(p, 3 * r, black_pen, steelblue_brush)
        D.Point(p, 3, red_pen, red_brush)

    # Flush save and show.
    D.FSS(filename = filename)

#---------------------------------------------------------------------------------------------------

def draw_ierarchical_tree(it,
                          deltas = (10, 40),
                          margins = (12, 12),
                          pen = aggdraw.Pen('orange', 2.0),
                          drawing_type = ClusteringDrawingType.Orthogonal,
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

        # Define line pen.
        line_pen = pen
        if it.IsOvershoot and ch.IsOvershoot:
            line_pen = aggdraw.Pen('black', 1.0)

        if drawing_type == ClusteringDrawingType.Lines:
            c.line(p + chp, line_pen)
        elif drawing_type == ClusteringDrawingType.Orthogonal:
            (px, py) = p
            (chpx, chpy) = chp
            c.line((px, py, chpx, py, chpx, chpy), line_pen)
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
    rp = lambda : (random.uniform(-50.0, 50.0), random.uniform(-50.0, 50.0))

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

def test_set_points2d_grid(nx, ny):
    """
    Test set with 2D points in grid nodes.

    Arguments:
        nx -- points number in X dimension,
        ny -- points number in Y dimension.

    Result:
        Points.
    """

    return fun.descartes_product(list(range(nx)), list(range(ny)))

#---------------------------------------------------------------------------------------------------

def test_set_points2d_circles(c, p):
    """
    Test set with 2D points in concentric circles.

    Arguments:
        c -- circles count,
        p -- points count.

    Result:
        Points.
    """

    ps = []

    for ci in range(c):
        r = 10.0 * ci
        da = 2.0 * math.pi / p
        ps = ps + [(r * math.cos(pi * da), r * math.sin(pi * da)) for pi in range(p)]

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

def test_case_points2d(ps, k, metric, metric_name, test_number = 1):
    """
    Run for test case.

    Arguments:
        ps -- points list,
        k -- clusters count,
        metric -- metric,
        metric_name -- name of metri,
        test_number -- number of test.
    """

    fmt = (metric_name, test_number)

    # Clustering and rearrange leafs Xs (because data is not ordered).
    tree = ierarchical_clustering(ps, k = k, metric = metric)
    tree.SetLeafsXs()
    tree.RefreshXs()

    # Find shoots_count.
    tree.CalculateHeightDifferences()
    tree.FindOvershoots(overshoots_count)

    # Ierarchical tree.
    #tree.Print()
    draw_ierarchical_tree(tree, filename = 'points2d_%s_tree_%d.png' % fmt)

    # Draw data.
    draw_data(tree,
              draw_clusters = False,
              pic_size = (600, 600), grid = (10.0, 10.0),
              filename = 'points2d_%s_init_%d.png' % fmt)
    draw_data(tree,
              draw_clusters = True,
              pic_size = (600, 600), grid = (10.0, 10.0),
              filename = 'points2d_%s_clusters_%d.png' % fmt)

    return tree.OvershootData()

#---------------------------------------------------------------------------------------------------

if __name__ == '__main__':

    clusters_count, overshoots_count = 12, 3
    run = RunType.Points2D

    if run == RunType.Points2D:

        # Test number.
        test_number = 1

        # Get test case.
        ps = test_set_points2d(200, k = clusters_count)
        #ps = test_set_points2d_grid(10, 10)
        #ps = test_set_points2d_circles(4, 50)

        overshoots = []
        modes = [(fun.partial_tail3(metric_tree_min_dist, 1.0), 'min1'),
                 (fun.partial_tail3(metric_tree_max_dist, 1.0), 'max1'),
                 (fun.partial_tail3(metric_tree_avg_dist, 1.0), 'avg1'),
                 (fun.partial_tail3(metric_tree_min_dist, 2.0), 'min2'),
                 (fun.partial_tail3(metric_tree_max_dist, 2.0), 'max2'),
                 (fun.partial_tail3(metric_tree_avg_dist, 2.0), 'avg2'),
                 (fun.partial_tail3(metric_tree_min_dist, 4.0), 'min4'),
                 (fun.partial_tail3(metric_tree_max_dist, 4.0), 'max4'),
                 (fun.partial_tail3(metric_tree_avg_dist, 4.0), 'avg4')]
        for (metric, metric_name) in modes:
            local_overshoots = test_case_points2d(ps, clusters_count,
                                                  metric = metric, metric_name = metric_name,
                                                  test_number = test_number)
            overshoots = overshoots + local_overshoots
            overshoots.sort()
            grouped_overshoots = lst.group(overshoots)
            sorted_overshoots = sorted([(c, p) for (p, c) in grouped_overshoots])
            sorted_overshoots.reverse()
            print("Grouped local overshoots for name %s:" % metric_name)
            print(sorted_overshoots)

        # Process statistics.
        draw_overshoots(ps,
                        sorted_overshoots,
                        overshoots_count,
                        pic_size = (600, 600),
                        grid = (10.0, 10.0),
                        filename = 'points2d_overshoots_%d.png' % test_number)

    elif run == RunType.Trajectory:
        ps = test_set_trajectory(0.03, 0.9)
        tree = ierarchical_clustering(ps, k = clusters_count)
        tree.CalculateHeightDifferences()
        tree.FindOvershoots(overshoots_count)
        tree.Print()
        draw_ierarchical_tree(tree,
                              filename = 'trajectory_tree.png')
        draw_data(tree,
                  draw_trajectory = True,
                  pic_size = (1000, 600),
                  grid = (1.0, 0.3),
                  filename = 'trajectory_init.png')
    else:
        raise Exception('wrong test run type')

#---------------------------------------------------------------------------------------------------
