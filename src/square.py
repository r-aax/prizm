# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 19:35:59 2019

@author: Rybakov
"""

import geom

a = geom.Vector(123.4, 348.6)
b = geom.Vector(84.4, 240.0)
c = geom.Vector(84.4, 233.6)
d = geom.Vector(92.0, 226.0)
e = geom.Vector(139.0, 208.0)
f = geom.Vector(144.4, 220.6)
g = geom.Vector(143.0, 228.6)
h = geom.Vector(143.0, 239.0)
i = geom.Vector(160.0, 318.0)
j = geom.Vector(160.0, 323.6)

lin = 57.5 / a.DistTo(b)

print("Points : ", [a, b, c, d, e, f, g, h, i, j])

aij = geom.Triangle(a, i, j).Area()
aib = geom.Triangle(a, i, b).Area()
ibh = geom.Triangle(i, b, h).Area()
bhc = geom.Triangle(b, h, c).Area()
chd = geom.Triangle(c, h, d).Area()
hdg = geom.Triangle(h, d, g).Area()
gdf = geom.Triangle(g, d, f).Area()
fde = geom.Triangle(f, d, e).Area()

sq = [aij, aib, ibh, bhc, chd, hdg, gdf, fde]

print("Relative squares : ", sq)

s = sum(sq) * lin * lin

print("Total square : ", s)

print(a.DistTo(j) * lin * 5.0 * 0.5)