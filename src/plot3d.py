# -*- coding: utf-8 -*-
"""
Plot3D supporting.

Created on Wed Jun  5 11:18:22 2019

@author: Rybakov
"""

#---------------------------------------------------------------------------------------------------
# Functions.
#---------------------------------------------------------------------------------------------------

def load_surface_peaces(filename):
    """
    Load peaces of surface.

    Arguments:
        filename -- name of file.

    Result:
        Loaded data.
    """

    # Open file.
    f = open(filename, 'r')

    # Surfaces peaces list.
    ss = []
    d1, d2, ps = 0, 0, []

    for l in f.readlines():

        if l == '\n':

            # If there is a peace of surface then add it to the list.
            if ps != []:
                ss.append((d1, d2, ps))
                d1, d2, ps = 0, 0, []

        elif l.startswith('PEACE_OF_SURFACE:'):

            # Init new peace of surface.
            t = l.split()
            d1, d2, ps = int(t[1][:-1]), int(t[2]), []

        else:

            # The line contains pooint.
            t = l.split()
            ps.append((float(t[0][1:-1]), float(t[1][:-1]), float(t[2][:-2])))

    f.close()

    return ss

#---------------------------------------------------------------------------------------------------

def save_plot3d(peaces, filename):
    """
    Save surface peaces to plot3d file.

    Arguments:
        peaces -- surface peaces,
        filename -- name of file.
    """

    # Count nodes.
    nodes_count = sum([len(ps) for (d1, d2, ps) in peaces])
    nodes_count2 = sum([d1 * d2 for (d1, d2, ps) in peaces])
    assert nodes_count == nodes_count2, 'wrong data in surface peaces file'

    # Count rectangles.
    triangles_count = 0
    for peace in peaces:
        (d1, d2, ps) = peace
        triangles_count += (d1 - 1) * (d2 - 1)
    triangles_count = triangles_count * 2

    f = open(filename, 'w')

    f.write('TITLE = "FE Surface Data ASCII"\n')
    f.write('VARIABLES = "X", "Y", "Z"\n')
    f.write('ZONE T="TRIANGLES", NODES="%d", ELEMENTS="%d", DATAPACKING="BLOCK", ZONETYPE="FETRIANGLE"\n'
            % (nodes_count, triangles_count))

    # Collect coordinates.
    xs = []
    ys = []
    zs = []
    for peace in peaces:
        (d1, d2, ps) = peace
        xs = xs + [x for (x, y, z) in ps]
        ys = ys + [y for (x, y, z) in ps]
        zs = zs + [z for (x, y, z) in ps]

    # Print coordinates.
    for x in xs:
        f.write('%f ' % x)
    f.write('\n')
    for y in ys:
        f.write('%f ' % y)
    f.write('\n')
    for z in zs:
        f.write('%f ' % z)
    f.write('\n')

    # Print links.
    off = 0
    for peace in peaces:
        (d1, d2, ps) = peace
        for i in range(d1 - 1):
            for j in range(d2 - 1):
                f.write('%d %d %d\n'
                        % ((off + 1) + i * d2 + j,
                           (off + 1) + (i + 1) * d2 + j,
                           (off + 1) + i * d2 + (j + 1)))
                f.write('%d %d %d\n'
                        % ((off + 1) + (i + 1) * d2 + j + 1,
                           (off + 1) + (i + 1) * d2 + j,
                           (off + 1) + i * d2 + (j + 1)))
        off = off + d1 * d2

    f.close()

#---------------------------------------------------------------------------------------------------
# Tests.
#---------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    peaces = load_surface_peaces('export1.txt')
    save_plot3d(peaces, 'air_inlet_01.dat')

#---------------------------------------------------------------------------------------------------
