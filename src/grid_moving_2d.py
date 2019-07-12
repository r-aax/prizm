# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 14:24:40 2019

@author: Rybakov
"""

import math
import geom
import draw
import aggdraw
import lst
import random
import vis

#---------------------------------------------------------------------------------------------------
# Functions.
#---------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------

def normal_1(v1, v2):
    """
    Simple realization of normal for this case.

    Arguments:
        v1 -- vector,
        v2 -- vector.

    Result:
        Normal.
    """

    o1 = v1.Orthogonal()
    if o1.Y < 0.0:
        o1.Negate()
    o2 = v2.Orthogonal()
    if o2.Y < 0.0:
        o2.Negate()
    o1.Normalize()
    o2.Normalize()
    n = o1 + o2
    n.Normalize()

    return n


def mod_dhs(a):
    """
    """

    #pos_a = sum([1 for ai in a if ai > 0])
    #neg_a = sum([1 for ai in a if ai < 0])

    #if pos_a > neg_a:
    #    return [max(ai, 0) for ai in a]
    #else:
    #    return [min(ai, 0) for ai in a]

    i = 0
    v = 0
    for j in range(1, len(a)):
        if abs(a[j]) > abs(v):
            v = a[j]
            i = j

    r = [0.0] * len(a)
    r[i] = v

    return r

#---------------------------------------------------------------------------------------------------
# Tests.
#---------------------------------------------------------------------------------------------------

def test(ss, case):

    # Generate points.
    # xs -- x coordinates of nodes
    # ps -- nodes points
    # cns -- cells normals
    # ns -- nodes normals
    # ss -- squares
    
    cells = len(ss)
    xs = [2.0 * math.pi * (i / cells) for i in range(cells + 1)]
    ps = [geom.Vector(x, math.sin(x)) for x in xs]
    cns = [None for _ in range(cells)]
    ns = [None for _ in range(cells + 1)]


    # Normals.
    for i, cn in enumerate(cns):
        cns[i] = normal_1(ps[i + 1] - ps[i], ps[i + 1] - ps[i])
    for i, n in enumerate(ns):
        if i == 0:
            ns[i] = normal_1(ps[i + 1] - ps[i], ps[i + 1] - ps[i])
        elif i == cells:
            ns[i] = normal_1(ps[i - 1] - ps[i], ps[i - 1] - ps[i])
        else:
            ns[i] = normal_1(ps[i - 1] - ps[i], ps[i + 1] - ps[i])

    D = draw.Drawer(draw_area = (-0.2 * math.pi, -1.1, 2.2 * math.pi, 2.1),
                    pic_size = (round(2.2 * math.pi * 100), round(3.2 * 100)))
    for i, p in enumerate(ps):
        if i < cells:
            D.Line(ps[i].Tuple(2), ps[i + 1].Tuple(2), pen = aggdraw.Pen('black', 3))
    for i, n in enumerate(ns):
        D.Line(ps[i].Tuple(2), (ps[i] + n).Tuple(2))
    for p in ps:
        D.Point(p.Tuple(2), 5, brush = aggdraw.Brush('black'))

    #
    #
    
    # Find alfa and beta.
    alfas = [None] * cells
    betas = [None] * cells

    ort_square = lambda i, h: h * ps[i].DistTo(ps[i + 1])
    square_i = lambda i, d1, d2: geom.Triangle(ps[i], ps[i + 1], ps[i] + ns[i].Scaled(d1)).Area() \
                             + geom.Triangle(ps[i + 1], ps[i] + ns[i].Scaled(d1), ps[i + 1] + ns[i + 1].Scaled(d2)).Area()
    square_h = lambda i, h: square_i(i, h / cns[i].CosAngle(ns[i]), h / cns[i].CosAngle(ns[i + 1]))


        #print(str(alfas))
        #print(str(betas))

    for i in range(cells):
        al = cns[i].CosAngle(ns[i])
        be = cns[i].CosAngle(ns[i + 1])
            #print(al, be)
        alc = math.acos(al)
        bec = math.acos(be)
            #print(alc, bec)
        if square_h(i, 0.001) > ort_square(i, 0.001):
            alfas[i] = 0.5 * math.pi + alc
            betas[i] = 0.5 * math.pi + bec
        else:
            alfas[i] = 0.5 * math.pi - alc
            betas[i] = 0.5 * math.pi - bec
    

    hs = [0.0] * len(ps)
    new_ps = [None for _ in ps]
    fun_l = lambda i: ps[i].DistTo(ps[i + 1])
    fun_t = lambda i: fun_l(i) * ss[i]
    fun_s = lambda i: 0.5 * (fun_l(i) * (hs[i] * math.sin(alfas[i]) + hs[i + 1] * math.sin(betas[i])) \
                             + hs[i] * hs[i + 1] * math.sin(alfas[i] + betas[i]))
    fun_delta = lambda i: math.pow(fun_s(i) - fun_t(i), 2.0)
    fun_diff = lambda i: 100.0 * (fun_s(i) - fun_t(i)) / fun_t(i)
    fun_full_t = lambda: sum([fun_l(i) * ss[i] for i in range(cells)])
    fun_full_s = lambda: sum([fun_s(i) for i in range(cells)])
    fun_full_delta = lambda: sum([fun_delta(i) for i in range(cells)])
    fun_full_dabs = lambda: sum([abs(fun_s(i) - fun_t(i)) for i in range(cells)])

    if case == 1:

        # PRISM case.

        for i in range(cells):
            if i == 0:
                d = cns[i].Scaled(2.0 * ss[i])
                D.Line(ps[i + 1].Tuple(2), (ps[i + 1] + d).Tuple(2))
                D.Line(ps[i].Tuple(2), (ps[i + 1] + d).Tuple(2))
            elif i == cells - 1:
                d = cns[i].Scaled(2.0 * ss[i])
                D.Line(ps[i].Tuple(2), (ps[i] + d).Tuple(2))
                D.Line((ps[i] + d).Tuple(2), ps[i + 1].Tuple(2))
            else:
                d = cns[i].Scaled(ss[i])
                D.Line(ps[i].Tuple(2), (ps[i] + d).Tuple(2))
                D.Line(ps[i + 1].Tuple(2), (ps[i + 1] + d).Tuple(2))
                D.Line((ps[i] + d).Tuple(2), (ps[i + 1] + d).Tuple(2))

        prev = ps[0]
        new_ps[0] = ps[0]
        for i in range(1, cells):
            if (i == 1) and (i == cells - 1):
                v = ss[i - 1] * 2.0
            elif i == 1:
                v = 0.5 * (2.0 * ss[i - 1] + ss[i])
            elif i == cells - 1:
                v = 0.5 * (ss[i - 1] + 2.0 * ss[i])
            else:
                v = 0.5 * (ss[i - 1] + ss[i])
            d = ns[i].Scaled(v)
            new_ps[i] = ps[i] + d
            D.Point(new_ps[i].Tuple(2), 3, brush = aggdraw.Brush('black'))
            D.Line(new_ps[i].Tuple(2), prev.Tuple(2), pen = aggdraw.Pen('black', 2))
            prev = new_ps[i]
        D.Line(prev.Tuple(2), ps[-1].Tuple(2), pen = aggdraw.Pen('black', 2))
        new_ps[-1] = ps[-1]

        hs = [ps[i].DistTo(new_ps[i]) for i in range(len(ps))]

        #print(hs)

        ress = (100.0 * (fun_full_s() - fun_full_t()) / fun_full_t(), 100.0 * fun_full_dabs() / fun_full_t())
        #print('RECT : ', ress)
        #return ress
        return [fun_diff(i) for i in range(cells)]

    elif case == 2:

        # PYRAMID case.

        true_square = lambda i: ss[i] * ps[i].DistTo(ps[i + 1])
        square_first = lambda d: geom.Triangle(ps[0], ps[1], ps[1] + ns[1].Scaled(d)).Area()
        square_last = lambda d: geom.Triangle(ps[-1], ps[-2], ps[-2] + ns[-2].Scaled(d)).Area()
        square_i = lambda i, d1, d2: geom.Triangle(ps[i], ps[i + 1], ps[i] + ns[i].Scaled(d1)).Area() \
                                     + geom.Triangle(ps[i + 1], ps[i] + ns[i].Scaled(d1), ps[i + 1] + ns[i + 1].Scaled(d2)).Area()
        square_h = lambda i, h: square_i(i, h / cns[i].CosAngle(ns[i]), h / cns[i].CosAngle(ns[i + 1]))


        lefts = [None for _ in range(cells)]
        rights = [None for _ in range(cells)]
        for i in range(cells):
            if i == 0:
                k = true_square(i) / square_first(0.001)
                d = ns[1].Scaled(0.001 * k)
                D.Line(ps[i + 1].Tuple(2), (ps[i + 1] + d).Tuple(2))
                D.Line(ps[i].Tuple(2), (ps[i + 1] + d).Tuple(2))
                rights[i] = ps[i + 1] + d
            elif i == cells - 1:
                k = true_square(i) / square_last(0.001)
                d = ns[-2].Scaled(0.001 * k)
                D.Line(ps[i].Tuple(2), (ps[i] + d).Tuple(2))
                D.Line(ps[i + 1].Tuple(2), (ps[i] + d).Tuple(2))
                lefts[i] = ps[i] + d
            else:
                s1 = square_h(i, 0.001)
                s2 = square_h(i, 0.002)
                b = (s1 - s2) / (0.001 - 0.002)
                a = s1 - b * 0.001
                h = (true_square(i) - a) / b
                d1 = ns[i].Scaled(h / cns[i].CosAngle(ns[i]))
                d2 = ns[i + 1].Scaled(h / cns[i].CosAngle(ns[i + 1]))
                D.Line((ps[i] + d1).Tuple(2), (ps[i + 1] + d2).Tuple(2))
                lefts[i] = ps[i] + d1
                rights[i] = ps[i + 1] + d2

        new_ps[0] = ps[0]
        new_ps[-1] = ps[-1]
        for i in range(1, cells):
            new_ps[i] = (rights[i - 1] + lefts[i]).Scaled(0.5)

        for i in range(len(new_ps) - 1):
            D.Line(new_ps[i].Tuple(2), new_ps[i + 1].Tuple(2), pen = aggdraw.Pen('black', 2))
        for p in new_ps:
            D.Point(p.Tuple(2), 3, brush = aggdraw.Brush('black'))

        hs = [ps[i].DistTo(new_ps[i]) for i in range(len(ps))]

        #print(hs)

        ress = (100.0 * (fun_full_s() - fun_full_t()) / fun_full_t(), 100.0 * fun_full_dabs() / fun_full_t())
        #print('TRAP : ', ress)
        #return ress
        return [fun_diff(i) for i in range(cells)]

    else:

        # Gradiaent descent.
        for i in range(len(ps)):
            new_ps[i] = ps[i]

        # ====================================
        fun_delta_der_next = lambda i: (fun_s(i) - fun_t(i)) * (fun_l(i) * math.sin(alfas[i]) - hs[i + 1] * math.sin(alfas[i] + betas[i]) )
        fun_delta_der_prev = lambda i: (fun_s(i - 1) - fun_t(i - 1)) * (fun_l(i - 1) * math.sin(betas[i - 1]) - hs[i - 1] * math.sin(alfas[i - 1] + betas[i - 1]))
        fun_delta_der = lambda i: (fun_delta_der_prev(i) + fun_delta_der_next(i))

        #print('before : ', fun_full_t(), fun_full_s())

        hs[0] = 0.0
        hs[-1] = 0.0
        for i in range(1, len(ps) - 1):
            if (i == 1) and (i == len(ps) - 2):
                hs[i] = 2.0 * ss[0]
            elif i == 1:
                hs[i] = 0.5 * (2.0 * ss[0] + ss[1])
            elif i == len(ps) - 2:
                hs[i] = 0.5 * (2.0 * ss[-1] + ss[-2])
            else:
                hs[i] = 0.5 * (ss[i - 1] + ss[i])

        perc_one = 0.02 * fun_full_t()
        #while fun_full_dabs() > perc_one:
        for _ in range(1000):
            dhs = [0.0] * len(ps)
            for i in range(1, len(ps) - 1):
                dhs[i] = fun_delta_der(i)
            #dhs = mod_dhs(dhs)
            #print(dhs)

            tau = 0.01
            for i in range(1, len(ps) - 1):
                ch = -tau * dhs[i]
                if (hs[i] + ch <= max(ss[i - 1], ss[i])) and (hs[i] + ch >= min(ss[i - 1], ss[i])):
                    hs[i] -= tau * dhs[i]
                else:
                    pass
                    #print('ref : ', hs[i], ch, ss[i - 1], ss[i])

            #print(' : iter : ', fun_full_t(), fun_full_s(), fun_full_dabs(), perc_one)

        #print(hs)

        ress = (100.0 * (fun_full_s() - fun_full_t()) / fun_full_t(), 100.0 * fun_full_dabs() / fun_full_t())
        #print('GRAD : ', ress)
        #return ress
        return [fun_diff(i) for i in range(cells)]


        # ====================================

        # ps
        for i in range(len(hs)):
            new_ps[i] = ps[i] + ns[i].Scaled(hs[i])

        for i in range(len(new_ps) - 1):
            D.Line(new_ps[i].Tuple(2), new_ps[i + 1].Tuple(2), pen = aggdraw.Pen('black', 2))
        for p in new_ps:
            D.Point(p.Tuple(2), 3, brush = aggdraw.Brush('black'))

    D.FSS('test.png')

#---------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    c = 1000
    dx = 2.0 * math.pi / c
    y = dx / 2.0
    #ss = [random.uniform(0.0, 2.0 * y) for _ in range(c)]
    ss = [y] * c
    ss[0] *= 0.5
    ss[-1] *= 0.5
    #print('ss = ', ss)
    r1 = test(ss, 1)
    r2 = test(ss, 2)
    vis.simple_graphic_ys(r1)
    vis.simple_graphic_ys(r2)
    
    print(r2)
    
    
