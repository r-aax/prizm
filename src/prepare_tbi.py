# -*- coding: utf-8 -*-
"""
Created on Sun Oct  6 11:47:04 2019

@author: Rybakov
"""

from functools import reduce

with open('../data/tri_box_intersect/tbi.txt', 'r') as f:

    ax, ay, az, bx, by, bz, cx, cy, cz, xl, xh, yl, yh, zl, zh, r = [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []

    l = f.readline()
    while l:

        ss = l.split()
        ax.append(ss[0])
        ay.append(ss[1])
        az.append(ss[2])
        bx.append(ss[3])
        by.append(ss[4])
        bz.append(ss[5])
        cx.append(ss[6])
        cy.append(ss[7])
        cz.append(ss[8])
        xl.append(ss[9])
        xh.append(ss[10])
        yl.append(ss[11])
        yh.append(ss[12])
        zl.append(ss[13])
        zh.append(ss[14])
        r.append(ss[15])

        l = f.readline()

    fax = open('../data/tri_box_intersect/ax.txt', 'w')
    fax.write(reduce(lambda x, y: x + ',' + y, ax))
    fax.close()

    fay = open('../data/tri_box_intersect/ay.txt', 'w')
    fay.write(reduce(lambda x, y: x + ',' + y, ay))
    fay.close()

    faz = open('../data/tri_box_intersect/az.txt', 'w')
    faz.write(reduce(lambda x, y: x + ',' + y, az))
    faz.close()

    fbx = open('../data/tri_box_intersect/bx.txt', 'w')
    fbx.write(reduce(lambda x, y: x + ',' + y, bx))
    fbx.close()

    fby = open('../data/tri_box_intersect/by.txt', 'w')
    fby.write(reduce(lambda x, y: x + ',' + y, by))
    fby.close()

    fbz = open('../data/tri_box_intersect/bz.txt', 'w')
    fbz.write(reduce(lambda x, y: x + ',' + y, bz))
    fbz.close()

    fcx = open('../data/tri_box_intersect/cx.txt', 'w')
    fcx.write(reduce(lambda x, y: x + ',' + y, cx))
    fcx.close()

    fcy = open('../data/tri_box_intersect/cy.txt', 'w')
    fcy.write(reduce(lambda x, y: x + ',' + y, cy))
    fcy.close()

    fcz = open('../data/tri_box_intersect/cz.txt', 'w')
    fcz.write(reduce(lambda x, y: x + ',' + y, cz))
    fcz.close()

    fxl = open('../data/tri_box_intersect/xl.txt', 'w')
    fxl.write(reduce(lambda x, y: x + ',' + y, xl))
    fxl.close()

    fxh = open('../data/tri_box_intersect/xh.txt', 'w')
    fxh.write(reduce(lambda x, y: x + ',' + y, xh))
    fxh.close()

    fyl = open('../data/tri_box_intersect/yl.txt', 'w')
    fyl.write(reduce(lambda x, y: x + ',' + y, yl))
    fyl.close()

    fyh = open('../data/tri_box_intersect/yh.txt', 'w')
    fyh.write(reduce(lambda x, y: x + ',' + y, yh))
    fyh.close()

    fzl = open('../data/tri_box_intersect/zl.txt', 'w')
    fzl.write(reduce(lambda x, y: x + ',' + y, zl))
    fzl.close()

    fzh = open('../data/tri_box_intersect/zh.txt', 'w')
    fzh.write(reduce(lambda x, y: x + ',' + y, zh))
    fzh.close()

    fr = open('../data/tri_box_intersect/r.txt', 'w')
    fr.write(reduce(lambda x, y: x + ',' + y, r))
    fr.close()

    f.close()