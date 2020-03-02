# -*- coding: utf-8 -*-
"""
Created on Sun Feb  9 15:55:41 2020

@author: Rybakov
"""

import vis
from draw import Drawer
import aggdraw
from PIL import Image, ImageDraw, ImageFont

#---------------------------------------------------------------------------------------------------
# Functions.
#---------------------------------------------------------------------------------------------------

def get_creatures_life_years():
    filename = '../../jdcsharp/NNBrothConsole/bin/Debug/result_01.txt'
    with open(filename, 'r') as f:
        li = []
        l = f.readline()
        while l:
            if 'Kill creature' in l:
                ss = l.split()
                li.append((int(ss[5][:-1]), int(ss[12][1:]), int(ss[14][:-1])))
            l = f.readline()
        return li

#---------------------------------------------------------------------------------------------------

def tr(s):
    while s[0] == '0':
        s = s[1:]
    return float(s)

#---------------------------------------------------------------------------------------------------

def extract_data():
    filename = '../../jdcsharp/NNBrothConsole/bin/Debug/result_01.txt'
    with open(filename, 'r') as f:
        li = []
        l = f.readline()
        while l:
            if 'Info' in l:
                ss = l.split()
                ss2 = ss[6][3:]
                ss3 = ss2.split('/')
                #
                ag = tr(ss[2][3:-1])
                lr = tr(ss[3][3:-1])
                rp = tr(ss[4][3:-1])
                tc = tr(ss[5][3:-1])
                nd = tr(ss3[0])
                ln = tr(ss3[1][:-2])
                #
                li.append((ag, lr, rp, tc, nd, ln))
            l = f.readline()
        return li

#---------------------------------------------------------------------------------------------------

def min_min_max_max_from_life_years(years):
    (i, b, d) = years[0]
    min_id, max_id, min_year, max_year = i, i, min(b, d), max(b, d)
    for (i, b, d) in years[1:]:
        min_id = min(min_id, i)
        max_id = max(max_id, i)
        min_year = min(min_year, min(b, d))
        max_year = max(max_year, max(b, d))
    return (min_id, max_id, min_year, max_year)

#---------------------------------------------------------------------------------------------------

def draw_life_years_diagram(years, bounds):
    max_len = 0
    for (i, b, d) in years:
        max_len = max(max_len, d - b + 1)
    min_id, max_id, min_year, max_year = bounds
    D = Drawer(draw_area = (0.0, 0.0, max_year, max_id),
               pic_size = (max_year, max_id),
               margins = (5, 5))
    for (i, b, d) in years:
        ll = d - b + 1
        if i == 1012:
            p = aggdraw.Pen('red', 10.0)
        elif ll > 100:
            p = aggdraw.Pen('green')
        elif ll > 30:
            p = aggdraw.Pen('lightblue')
        else:
            p = aggdraw.Pen('gray')
        print(max_len)
        D.Line((b, i), (d, i), pen = p)
    D.Rect((0, 0), (max_year, max_id), pen = p)
    D.FSS()

#---------------------------------------------------------------------------------------------------
# Tests.
#---------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    #
    #life_years = get_creatures_life_years()
    #bounds = min_min_max_max_from_life_years(life_years)
    #print(bounds)
    #draw_life_years_diagram(life_years, bounds)
    #
    data = extract_data()
    ags = [ag for (ag, lr, rp, tc, nd, ln) in data]
    lrs = [lr for (ag, lr, rp, tc, nd, ln) in data]
    rps = [rp for (ag, lr, rp, tc, nd, ln) in data]
    tcs = [tc for (ag, lr, rp, tc, nd, ln) in data]
    nds = [nd for (ag, lr, rp, tc, nd, ln) in data]
    lns = [ln for (ag, lr, rp, tc, nd, ln) in data]
    nn = 6
    if nn == 1:
        vis.simple_graphic_ys(ags,
                              title = 'Возраст лидера популяции',
                              x_label = 'номер поколения',
                              y_label = 'возраст')
    elif nn == 2:
        vis.simple_graphic_ys(lrs,
                              title = 'Количество итераций обучения лидера популяции',
                              x_label = 'номер поколения',
                              y_label = 'количество итераций обучения')
    elif nn == 3:
        vis.simple_graphic_ys(rps,
                              title = 'Доля верных ответов лидера популяции',
                              x_label = 'номер поколения',
                              y_label = 'доля верных ответов')
    elif nn == 4:
        vis.simple_graphic_ys(tcs,
                              title = 'Значение функции стоимости лидера популяции',
                              x_label = 'номер поколения',
                              y_label = 'значение функции стоимости')
    elif nn == 5:
        vis.simple_graphic_ys(nds,
                              title = 'Порядок графа лидера популяции',
                              x_label = 'номер поколения',
                              y_label = 'порядок')
    elif nn == 6:
        vis.simple_graphic_ys(lns,
                              title = 'Размер графа лидера популяции',
                              x_label = 'номер поколения',
                              y_label = 'размер')

#---------------------------------------------------------------------------------------------------
