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
import aggdraw
import mth
import metrics
from functools import reduce
from draw import Drawer
from enum import Enum
from PIL import Image
from htree import HTree

#---------------------------------------------------------------------------------------------------

def NewHTree(data = None):
    """
    Constructor.

    Arguments:
        data -- data.
    """

    # Create tree.
    ht = HTree()

    # X coordinate for visualization.
    ht.X = None

    # Mark for additional purposes.
    ht.Mark = False

    # Cluster number.
    ht.KN = -1

    # Data.
    ht.Data = data

    # Outlier value.
    ht.OutlierValue = None

    # Is outlier.
    ht.IsOutlier = False

    return ht

#---------------------------------------------------------------------------------------------------

def NextLeafLeftRoRight(ht, leaf):
    """
    Get next leaf (left to right walk).

    Arguments:
        ht -- hierarchical tree,
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

def ClusterNumber(ht):
    """
    Node cluster number.

    Arguments:
        ht -- hierarchical tree.

    Result:
        Cluster number.
    """

    kn = ht.KN

    if ht.IsRoot():
        return kn
    elif kn != -1:
        return kn
    else:
        return ClusterNumber(ht.Parent)

#---------------------------------------------------------------------------------------------------

def Str(ht):
    """
    String representation.

    Arguments:
        ht -- hierarchical tree.

    Result:
        String.
    """

    ov_str = (', hd = ' + str(ht.OutlierValue)) if ht.OutlierValue != None else ''

    if ht.Data == None:
        return 'None'
    elif isinstance(ht.Data, numbers.Number):
        return 'x = %f, kn = %d, data = %f%s (1d)' \
               % (ht.X, ht.KN, ht.Data, ov_str)
    elif isinstance(ht.Data, tuple):
        if len(ht.Data) == 2:
            return 'x = %f, kn = %d, data = (%f, %f)%s (2d)' \
                   % (ht.X, ht.KN, ht.Data[0], ht.Data[1], ov_str)
        else:
            raise Exception('wrong data for string representation')
    else:
        return 'str'

#---------------------------------------------------------------------------------------------------

def Print(ht):
    """
    Print tree.

    Arguments:
        ht -- hierarchical tree.
    """

    level = ht.Level()
    indent = ' ' * level
    print('%s [L%d] : %s' % (indent, level, Str(ht)))

    # Print all children.
    for c in ht.Children:
        Print(c)

#---------------------------------------------------------------------------------------------------

def TreeSizes(ht, deltas, margins):
    """
    Tree sizes.

    Arguments:
        ht -- hierarchical tree,
        deltas -- distances between nodes,
        margins -- margins.

    Result:
        Tree sizes.
    """

    dx, dy = deltas
    mx, my = margins

    return (int((ht.Width() - 1) * dx + 2 * mx), int((ht.Height() - 1) * dy + 2 * my))

#---------------------------------------------------------------------------------------------------

def NodeCoordinates(ht, deltas, margins):
    """
    Coordinates of node.

    Arguments:
        ht -- hierarchical tree,
        deltas -- distances between nodes,
        margins -- margins.

    Result:
        Node coordinates.
    """

    dx, dy = deltas
    mx, my = margins

    return (int(ht.X * dx + mx), int((ht.Height() - 1) * dy + my))

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
    t = NewHTree()

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

def CalculateOutlierValues(ht):
    """
    Calculate outlier values.

    Arguments:
        ht -- hierarchical tree.
    """

    # Calculate own heights difference.
    if ht.IsRoot():
        ht.OutlierValue = 0.0
    elif ClusterNumber(ht) == -1:
        ht.OutlierValue = 0.0
    else:
        ht.OutlierValue = float(ht.Parent.Height() - ht.Height()) / ht.Width()

    # Calculate childrens' heights differences.
    if not ht.IsLeaf():
        for ch in ht.Children:
            CalculateOutlierValues(ch)

#---------------------------------------------------------------------------------------------------

def MaxOutlierLeaf(ht):
    """
    Maximum outlier leaf.

    Arguments:
        ht -- hierarchical tree.

    Result:
        Maximum outlier leaf.
    """

    if ht.IsOutlier:
        return None
    elif ht.IsLeaf():
        return ht
    else:
        ms = [MaxOutlierLeaf(ch) for ch in ht.Children]
        cur_m = ht
        for m in ms:
           if cur_m == None:
               cur_m = m
           elif m == None:
               pass
           elif m.OutlierValue > cur_m.OutlierValue:
                cur_m = m
           else:
               pass
        return cur_m

#---------------------------------------------------------------------------------------------------

def SetOutlierRecursive(ht):
    """"
    Set outlier recursive.

    Arguments:
        ht -- hierarchical tree.
    """

    ht.IsOutlier = True

    for ch in ht.Children:
        SetOutlierRecursive(ch)

#---------------------------------------------------------------------------------------------------

def FindOutliers(ht, n):
    """
    Find outliers.

    Arguments:
        ht -- hierarchical tree,
        n -- count.
    """

    for i in range(n):
        m = MaxOutlierLeaf(ht)
        SetOutlierRecursive(m)

#---------------------------------------------------------------------------------------------------

def SetLeafsXs(ht, start = 0):
    """
    Set X coordinates for leafs.

    Arguments:
        ht -- hierarchical tree,
        start -- start position.
    """

    if ht.IsLeaf():
        ht.X = start
    else:
        cur = start
        for i in range(ht.ChildrenCount()):
            ch = ht.Children[i]
            SetLeafsXs(ch, cur)
            cur = cur + ch.Width()

#---------------------------------------------------------------------------------------------------

def RefreshXs(ht):
    """
    Refresh Xs.

    Arguments:
        ht -- hierarchical tree.
    """

    if ht.IsLeaf():
        pass
    else:
        for ch in ht.Children:
            RefreshXs(ch)
        xws = [(ch.X, ch.Width()) for ch in ht.Children]
        (xs, ws) = fun.unzip(xws)
        ht.X = mth.avg_weighted(xs, ws)

#---------------------------------------------------------------------------------------------------

def AllData(ht):
    """
    All data of tree (leafs).

    Arguments:
        ht -- hierarchical tree.

    Result:
        Data as a list.
    """

    if ht.IsLeaf():
        return [ht.Data]
    else:
        return reduce(operator.__concat__,
                      [AllData(ch) for ch in ht.Children],
                      [])

#---------------------------------------------------------------------------------------------------

def OutlierData(ht):
    """
    Outlierdata (leafs).

    Arguments:
        ht -- hierarchical tree.

    Result:
        Outlier data.
    """

    if ht.IsLeaf():
        if ht.IsOutlier:
            return [ht.Data]
        else:
            return []
    else:
        return reduce(operator.__concat__,
                      [OutlierData(ch) for ch in ht.Children],
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

def lp_norm_list(t1, t2, r = 2.0):
    """
    List of distances for all pair of leafs.

    Arguments:
        t1 -- first tree,
        t2 -- second tree,
        r -- power.

    Result:
        List of distances between nodes from two trees.
    """

    d1 = AllData(t1)
    d2 = AllData(t2)
    datas = lst.descartes_product(d1, d2)

    return [metrics.lp_norm(p1, p2, r) for (p1, p2) in datas]

#---------------------------------------------------------------------------------------------------

def sup_norm_list(t1, t2):
    """
    List of supreme norm of leafs.

    Arguments:
        t1 -- first tree,
        t2 -- second tree.

    Result:
        List of supreme norms.
    """

    d1 = AllData(t1)
    d2 = AllData(t2)
    datas = lst.descartes_product(d1, d2)

    return [metrics.sup_norm(p1, p2) for (p1, p2) in datas]

#---------------------------------------------------------------------------------------------------

def jeffreys_matsushita_list(t1, t2):
    """
    List of Jeffreys-Matsushita.

    Arguments:
        t1 -- first tree,
        t2 -- second tree.

    Result:
        List of Jeffreys-Matsushita.
    """

    d1 = AllData(t1)
    d2 = AllData(t2)
    datas = lst.descartes_product(d1, d2)

    return [metrics.jeffreys_matsushita(p1, p2) for (p1, p2) in datas]

#---------------------------------------------------------------------------------------------------

def div_coef_list(t1, t2):
    """
    List of divergences coefficients.

    Arguments:
        t1 -- first tree,
        t2 -- second tree.

    Result:
        List of divergences coefficients.
    """

    d1 = AllData(t1)
    d2 = AllData(t2)
    datas = lst.descartes_product(d1, d2)

    return [metrics.div_coef(p1, p2) for (p1, p2) in datas]

#---------------------------------------------------------------------------------------------------

def metric_tree_min_lp_norm(t1, t2, r = 2.0):
    """
    Metric - minimal distance between nodes.

    Arguments:
        t1 -- first tree,
        t2 -- second tree,
        r -- power.

    Result:
        Metric result.
    """

    return min(lp_norm_list(t1, t2, r))

#---------------------------------------------------------------------------------------------------

def metric_tree_max_lp_norm(t1, t2, r = 2.0):
    """
    Metric - maximal distance between nodes.

    Arguments:
        t1 -- first tree,
        t2 -- second tree,
        r -- power.

    Result:
        Metric result.
    """

    return max(lp_norm_list(t1, t2, r))

#---------------------------------------------------------------------------------------------------

def metric_tree_avg_lp_norm(t1, t2, r = 2.0):
    """
    Metric - average distance between nodes.

    Arguments:
        t1 -- first tree,
        t2 -- second tree,
        r -- power.

    Result:
        Metric result.
    """

    return mth.avg_arith(lp_norm_list(t1, t2, r))

#---------------------------------------------------------------------------------------------------

def metric_tree_min_sup_norm(t1, t2):
    """
    Metric - minimal supreme norm.

    Arguments:
        t1 -- first tree,
        t2 -- second tree.

    Result:
        Metric result.
    """

    return min(sup_norm_list(t1, t2))

#---------------------------------------------------------------------------------------------------

def metric_tree_max_sup_norm(t1, t2):
    """
    Metric - maximal supreme norm.

    Arguments:
        t1 -- first tree,
        t2 -- second tree.

    Result:
        Metric result.
    """

    return max(sup_norm_list(t1, t2))

#---------------------------------------------------------------------------------------------------

def metric_tree_avg_sup_norm(t1, t2):
    """
    Metric - average supreme norm.

    Arguments:
        t1 -- first tree,
        t2 -- second tree.

    Result:
        Metric result.
    """

    return mth.avg_arith(sup_norm_list(t1, t2))

#---------------------------------------------------------------------------------------------------

def metric_tree_min_jeffreys_matsushita(t1, t2):
    """
    Metric - minimal Jeffreys Matsushita between nodes.

    Arguments:
        t1 -- first tree,
        t2 -- second tree.

    Result:
        Metric result.
    """

    return min(jeffreys_matsushita_list(t1, t2))

#---------------------------------------------------------------------------------------------------

def metric_tree_max_jeffreys_matsushita(t1, t2):
    """
    Metric - maximal Jeffreys Matsushita between nodes.

    Arguments:
        t1 -- first tree,
        t2 -- second tree.

    Result:
        Metric result.
    """

    return max(jeffreys_matsushita_list(t1, t2))

#---------------------------------------------------------------------------------------------------

def metric_tree_avg_jeffreys_matsushita(t1, t2):
    """
    Metric - average Jeffreys Matsushita between nodes.

    Arguments:
        t1 -- first tree,
        t2 -- second tree.

    Result:
        Metric result.
    """

    return mth.avg_arith(jeffreys_matsushita_list(t1, t2))

#---------------------------------------------------------------------------------------------------

def metric_tree_min_div_coef(t1, t2):
    """
    Metric - minimal divergence coefficient between nodes.

    Arguments:
        t1 -- first tree,
        t2 -- second tree.

    Result:
        Metric result.
    """

    return min(div_coef_list(t1, t2))

#---------------------------------------------------------------------------------------------------

def metric_tree_max_div_coef(t1, t2):
    """
    Metric - maximal divergence coefficient between nodes.

    Arguments:
        t1 -- first tree,
        t2 -- second tree.

    Result:
        Metric result.
    """

    return max(div_coef_list(t1, t2))

#---------------------------------------------------------------------------------------------------

def metric_tree_avg_div_coef(t1, t2):
    """
    Metric - average divergence coefficient between nodes.

    Arguments:
        t1 -- first tree,
        t2 -- second tree.

    Result:
        Metric result.
    """

    return mth.avg_arith(div_coef_list(t1, t2))

#---------------------------------------------------------------------------------------------------

def hierarchical_clustering(ps,
                            k = 1,
                            metric = metric_tree_avg_lp_norm,
                            nearest_type = ClusteringNearestType.All):
    """
    Hierarchical clustering.

    Arguments:
        ps -- array of points,
        k -- clusters count for mark,
        metric -- metric function.

    Result:
        Clustering tree.
    """

    trees = [NewHTree(p) for p in ps]

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
        new_tree = Merge(trees[fpi], trees[spi])
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
            leaf2 = NextLeafLeftRoRight(tree, leaf1)
            while leaf2 != None:

                # Draw.
                d1 = leaf1.Data
                d2 = leaf2.Data
                kn1 = ClusterNumber(leaf1)
                kn2 = ClusterNumber(leaf2)
                is_not_outliers = (not leaf1.IsOutlier) and (not leaf2.IsOutlier)
                if is_not_outliers and (kn1 != -1) and (kn1 == kn2):
                    D.Line(d1, d2, pen = aggdraw.Pen(pretty_color(kn1), 1.0))

                leaf2 = NextLeafLeftRoRight(tree, leaf2)
            leaf1 = NextLeafLeftRoRight(tree, leaf1)

    backcolor = D.Backcolor
    backcolor_pen = aggdraw.Pen(backcolor, 1.0)
    backcolor_brush = aggdraw.Brush(backcolor)

    # Draw points.
    leaf = tree.LeftLeaf()
    while leaf != None:

        # Default colors and etc.
        color = 'red'
        point_radius = 3

        # Define color if cluster number is set.
        if draw_clusters:
            kn = ClusterNumber(leaf)
            if leaf.IsOutlier:
                color = 'black'
                point_radius = 5
            elif kn != -1:
                color = pretty_color(kn)

        # Define pen, brush and draw point.
        point_pen = aggdraw.Pen(color, 1.0)
        point_brush = aggdraw.Brush(color)
        D.Point(leaf.Data, point_radius, pen = point_pen, brush = point_brush)
        if draw_clusters:
            if not leaf.IsOutlier:
                D.Point(leaf.Data, 1, pen = backcolor_pen, brush = backcolor_brush)
        leaf = NextLeafLeftRoRight(tree, leaf)

    # Flush save and show.
    D.FSS(filename = filename)

#---------------------------------------------------------------------------------------------------

def draw_outliers(ps, outliers, ok,
                  pic_size = (640, 480),
                  is_axis = True,
                  grid = None,
                  filename = None):
    """
    Draw outliers.

    Arguments:
        ps -- points,
        outliers -- list of outliers,
        ok -- outliers count (leaders),
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

    # Draw outliers.
    black_pen = aggdraw.Pen('black', 2.0)
    steelblue_brush = aggdraw.Brush('steelblue')
    for outlier in outliers[ok:]:
        (r, p) = outlier
        D.Point(p, 3 * r, black_pen)
    for outlier in outliers[:ok]:
        (r, p) = outlier
        D.Point(p, 3 * r, black_pen, steelblue_brush)
        D.Point(p, 3, red_pen, red_brush)

    # Flush save and show.
    D.FSS(filename = filename)

#---------------------------------------------------------------------------------------------------

def draw_hierarchical_tree(ht,
                           deltas = (10, 40),
                           margins = (12, 12),
                           pen = aggdraw.Pen('orange', 2.0),
                           drawing_type = ClusteringDrawingType.Orthogonal,
                           filename = None):
    """
    Draw hierarchical tree.

    Arguments:
        ht -- hierarchical tree,
        deltas -- distances between nodes,
        margins -- margins,
        pen -- pen,
        drawing_type -- drawing type.
    """

    # Create image.
    img = Image.new('RGB', TreeSizes(ht, deltas, margins), color = (255, 255, 255))
    c = aggdraw.Draw(img)
    c.setantialias(True)

    # Recursive draw.
    draw_hierarchical_tree_on_img(ht, c, deltas, margins, pen, drawing_type)

    # Flush, save and show.
    c.flush()
    if filename != None:
        img.save(filename)
    img.show()

#---------------------------------------------------------------------------------------------------

def draw_hierarchical_tree_on_img(ht, c, deltas, margins, pen, drawing_type):
    """
    Draw hierarchical tree on image.

    Arguments:
        ht -- hierarchical tree,
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
    if ht.Mark:
        pen = aggdraw.Pen(pretty_color(ht.KN), 2.0)

    # Coordinates.
    p = NodeCoordinates(ht, deltas, margins)

    # Draw children.
    for ch in ht.Children:
        chp = NodeCoordinates(ch, deltas, margins)

        # Define line pen.
        line_pen = pen
        if ht.IsOutlier and ch.IsOutlier:
            line_pen = aggdraw.Pen('black', 1.0)

        if drawing_type == ClusteringDrawingType.Lines:
            c.line(p + chp, line_pen)
        elif drawing_type == ClusteringDrawingType.Orthogonal:
            (px, py) = p
            (chpx, chpy) = chp
            c.line((px, py, chpx, py, chpx, chpy), line_pen)
        else:
            raise Exception('wrong drawing type : %s' % str(drawing_type))

        draw_hierarchical_tree_on_img(ch, c, deltas, margins, pen, drawing_type)

    # Draw point.
    if ht.IsOutlier:
        c.ellipse(expand_to_circle(p, 5), aggdraw.Pen('black', 2.0), aggdraw.Brush('black'))
    else:
        c.ellipse(expand_to_circle(p, 3), pen, brush)
    if ht.Mark:
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

def test_set_trajectory(ripple, outlier):
    """
    Test set Trajectory.

    Arguments:
        ripple -- ripple intensive,
        outlier -- outlier value.

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

    # Outliers.
    for i in range(4):
        ind = 40 * (i + 1)
        dx = random.uniform(-outlier, outlier)
        dy = random.uniform(-outlier, outlier)
        ps[ind] = (ps[ind][0] + dx, ps[ind][1] + dy)

    return ps

#---------------------------------------------------------------------------------------------------
# Run type.
#---------------------------------------------------------------------------------------------------

class RunType(Enum):
    """
    Type of test run.
    """

    # Simple test.
    Test = 0

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
    tree = hierarchical_clustering(ps, k = k, metric = metric)
    SetLeafsXs(tree)
    RefreshXs(tree)

    # Find shoots_count.
    CalculateOutlierValues(tree)
    FindOutliers(tree, outliers_count)

    # Ierarchical tree.
    #tree.Print()
    draw_hierarchical_tree(tree, filename = 'points2d_%s_tree_%d.png' % fmt)

    # Draw data.
    draw_data(tree,
              draw_clusters = False,
              pic_size = (600, 600), grid = (10.0, 10.0),
              filename = 'points2d_%s_init_%d.png' % fmt)
    draw_data(tree,
              draw_clusters = True,
              pic_size = (600, 600), grid = (10.0, 10.0),
              filename = 'points2d_%s_clusters_%d.png' % fmt)

    return OutlierData(tree)

#---------------------------------------------------------------------------------------------------

if __name__ == '__main__':

    points_count, clusters_count, outliers_count = 50, 6, 3
    run = RunType.Points2D

    if run == RunType.Test:

        # Simple test.
        pass

    elif run == RunType.Points2D:

        # Test number.
        test_number = 1

        # Get test case.
        ps = test_set_points2d(points_count, k = clusters_count)

        outliers = []
        modes = [
                 #(fun.partial_tail3(metric_tree_min_lp_norm, 1.0), 'min1'),
                 #(fun.partial_tail3(metric_tree_max_lp_norm, 1.0), 'max1'),
                 #(fun.partial_tail3(metric_tree_avg_lp_norm, 1.0), 'avg1'),
                 #(fun.partial_tail3(metric_tree_min_lp_norm, 2.0), 'min2'),
                 #(fun.partial_tail3(metric_tree_max_lp_norm, 2.0), 'max2'),
                 #(fun.partial_tail3(metric_tree_avg_lp_norm, 2.0), 'avg2'),
                 #(fun.partial_tail3(metric_tree_min_lp_norm, 4.0), 'min4'),
                 #(fun.partial_tail3(metric_tree_max_lp_norm, 4.0), 'max4'),
                 #(fun.partial_tail3(metric_tree_avg_lp_norm, 4.0), 'avg4'),
                 #(metric_tree_min_sup_norm, 'min_s'),
                 #(metric_tree_max_sup_norm, 'max_s'),
                 #(metric_tree_avg_sup_norm, 'avg_s'),
                 #(metric_tree_min_jeffreys_matsushita, 'min_jm'),
                 #(metric_tree_max_jeffreys_matsushita, 'max_jm'),
                 #(metric_tree_avg_jeffreys_matsushita, 'avg_jm'),
                 #(metric_tree_min_div_coef, 'min_dc'),
                 #(metric_tree_max_div_coef, 'max_dc'),
                 (metric_tree_avg_div_coef, 'avg_dc')
                 ]
        for (metric, metric_name) in modes:
            local_outliers = test_case_points2d(ps, clusters_count,
                                                metric = metric, metric_name = metric_name,
                                                test_number = test_number)
            outliers = outliers + local_outliers
            outliers.sort()
            grouped_outliers = lst.group(outliers)
            sorted_outliers = sorted([(c, p) for (p, c) in grouped_outliers])
            sorted_outliers.reverse()
            print('Sorted local outliers for name %s:' % metric_name)
            print(sorted_outliers)

        # Process statistics.
        sorted_outliers = [(c, p) for (c, p) in sorted_outliers if c > 1]
        print('Final sorted outliers:')
        print(sorted_outliers)
        draw_outliers(ps,
                        sorted_outliers,
                        outliers_count,
                        pic_size = (600, 600),
                        grid = (10.0, 10.0),
                        filename = 'points2d_outliers_%d.png' % test_number)

    elif run == RunType.Trajectory:

        # Not implemented yet.
        pass

    else:
        raise Exception('wrong test run type')

#---------------------------------------------------------------------------------------------------
