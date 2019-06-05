# -*- coding: utf-8 -*-
"""
Plot3D supporting.

Created on Wed Jun  5 11:18:22 2019

@author: Rybakov
"""

from gsu import Grid

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

    g = Grid()

    # Construct grid.
    for peace in peaces:
        (d1, d2, ps) = peace
        g.ConstructFromVectorsFlatMatrix(d1, d2, ps)

    # Export.
    g.ExportToTecplot(filename)

#---------------------------------------------------------------------------------------------------
# Tests.
#---------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    peaces = load_surface_peaces('faces_from_gridmaster.txt')
    save_plot3d(peaces, 'air_inlet_3.dat')

#---------------------------------------------------------------------------------------------------
