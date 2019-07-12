# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 15:18:49 2019

@author: Rybakov
"""

from functools import reduce
import lst
import gsu

#---------------------------------------------------------------------------------------------------
# Functions.
#---------------------------------------------------------------------------------------------------

def extract_data_from_dat(filename, ofilename, mask):
    """
    Extract data from *.dat file.

    Arguments:
        filename -- name of input file,
        ofilename -- name of output file,
        mask -- mask.
    """

    is_title_str = lambda s: 'TITLE' in s
    is_variables_str = lambda s: 'VARIABLES' in s
    is_zone_str = lambda s: 'ZONE' in s
    is_meta_str = lambda s: is_title_str(s) or is_variables_str(s) or is_zone_str(s)

    with open(filename, 'r') as f:
        with open(ofilename, 'w') as of:

            wait_nodes, wait_faces = 0, 0
            l = f.readline()
            while l:

                if is_meta_str(l):
                    if (wait_nodes > 0) or (wait_faces > 0):
                        raise Exception('internal error')
                    if is_title_str(l):
                        print('ORIGINAL : ', l)
                    elif is_variables_str(l):
                        print('ORIGINAL : ', l)
                        variables = l[l.index('=') + 1 : ].split(',')
                        variables_count = len(variables)
                        print('variables = ', variables, ' count = ', variables_count)
                    elif is_zone_str(l):
                        print('ORIGINAL : ', l)
                        fields = l[l.index('N=') : ].split()
                        nodes_count, faces_count = int(fields[1]), int(fields[4])
                        print('nodes_count = ', nodes_count, ' faces_count = ', faces_count)
                        wait_nodes, wait_faces = nodes_count, faces_count
                else:
                    if wait_nodes > 0:
                        ms = lst.mask(l.split(), mask)
                        if ms != []:
                            out = reduce(lambda a, b: a + ' ' + b, ms) + '\n'
                            of.write(out)
                        wait_nodes -= 1
                    elif wait_faces > 0:
                        wait_faces -= 1

                l = f.readline()

            of.close()
        f.close()

#---------------------------------------------------------------------------------------------------

def merge_files(f1name, f2name, outname):
    """
    Merge two files.

    Arguments:
        f1name -- first file,
        f2name -- second file,
        outname -- out file.
    """

    with open(f1name, 'r') as f1:
        with open(f2name, 'r') as f2:
            with open(outname, 'w') as out:

                l1 = f1.readline()
                l2 = f2.readline()

                while l1:

                    s1 = l1.split()
                    s2 = l2.split()[3:]
                    s = s1 + s2

                    v1 = float(s[3])
                    v2 = float(s[4])
                    v3 = float(s[5])
                    v4 = float(s[6])
                    v5 = float(s[7])

                    if (v1 != 0.0) or (v2 != 0.0) or (v3 != 0.0) or (v4 != 0.0) or (v5 != 0.0):
                        lout = reduce(lambda a, b: a + ' ' + b, s) + '\n'
                        out.write(lout)

                    l1 = f1.readline()
                    l2 = f2.readline()

                out.close()
            f2.close()
        f1.close()

#---------------------------------------------------------------------------------------------------

def step_0_cut_extra_variables(filename, ofilename):
    """
    Cut extra variables.

    Arguments:
        filename -- name of input file,
        ofilename -- name of output file.
    """

    is_title_str = lambda s: 'TITLE' in s
    is_variables_str = lambda s: 'VARIABLES' in s
    is_zone_str = lambda s: 'ZONE' in s
    is_meta_str = lambda s: is_title_str(s) or is_variables_str(s) or is_zone_str(s)

    with open(filename, 'r') as f:
        with open(ofilename, 'w') as of:

            # Write head to out file.
            of.write('# EXPORT MODE: CHECK_POINT\n')
            of.write('TITLE = "FE Surface Data ASCII"\n')
            of.write('VARIABLES = "X", "Y", "Z"\n')

            wait_nodes, wait_faces = 0, 0
            l = f.readline()
            while l:

                if is_meta_str(l):
                    if (wait_nodes > 0) or (wait_faces > 0):
                        raise Exception('internal error')
                    if is_title_str(l):
                        print('ORIGINAL : ', l)
                    elif is_variables_str(l):
                        print('ORIGINAL : ', l)
                        variables = l[l.index('=') + 1 : ].split(',')
                        variables_count = len(variables)
                        print('variables = ', variables, ' count = ', variables_count)
                    elif is_zone_str(l):
                        print('ORIGINAL : ', l)
                        fields = l[l.index('N=') : ].split()
                        nodes_count, faces_count = int(fields[1]), int(fields[4])
                        print('nodes_count = ', nodes_count, ' faces_count = ', faces_count)
                        wait_nodes, wait_faces = nodes_count, faces_count
                        of.write('ZONE T="TRIANGLES", NODES=%d, ELEMENTS=%d, DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL\n'
                                 % (nodes_count, faces_count))
                else:
                    if wait_nodes > 0:
                        ms = l.split()[:3]
                        of.write(reduce(lambda a, b: a + ' ' + b, ms) + '\n')
                        wait_nodes -= 1
                    elif wait_faces > 0:
                        of.write(l)
                        wait_faces -= 1

                l = f.readline()

            of.close()
        f.close()

#---------------------------------------------------------------------------------------------------

def step_0_pack_point_to_block(filename, ofilename):
    """
    Change data packing type.

    Arguments:
        filename -- input file name,
        ofilename -- output file name.
    """

    is_title_str = lambda s: 'TITLE' in s
    is_variables_str = lambda s: 'VARIABLES' in s
    is_zone_str = lambda s: 'ZONE' in s
    is_meta_str = lambda s: is_title_str(s) or is_variables_str(s) or is_zone_str(s)

    with open(filename, 'r') as f:
        with open(ofilename, 'w') as of:

            # Write head to out file.
            of.write('# EXPORT MODE: CHECK_POINT\n')
            of.write('TITLE = "FE Surface Data ASCII"\n')
            of.write('VARIABLES = "X", "Y", "Z"\n')

            wait_nodes, wait_faces = 0, 0
            storage = []
            l = f.readline()
            while l:

                if is_meta_str(l):
                    if (wait_nodes > 0) or (wait_faces > 0):
                        raise Exception('internal error')
                    if is_title_str(l):
                        print('ORIGINAL : ', l)
                    elif is_variables_str(l):
                        print('ORIGINAL : ', l)
                        variables = l[l.index('=') + 1 : ].split(',')
                        variables_count = len(variables)
                        print('variables = ', variables, ' count = ', variables_count)
                    elif is_zone_str(l):
                        print('ORIGINAL : ', l)
                        fields = l[l.index('NODES=') : ].split()
                        print('fields = ', fields)
                        nodes_count, faces_count = int(fields[0][6:-1]), int(fields[1][9:-1])
                        print('nodes_count = ', nodes_count, ' faces_count = ', faces_count)
                        wait_nodes, wait_faces = nodes_count, faces_count
                        of.write('ZONE T="TRIANGLES", NODES=%d, ELEMENTS=%d, DATAPACKING=BLOCK, ZONETYPE=FETRIANGLE\n'
                                 % (nodes_count, faces_count * 2))
                else:
                    if wait_nodes > 0:
                        ms = l.split()[:3]
                        storage.append(ms)
                        wait_nodes -= 1

                        if wait_nodes == 0:
                            for [a, _, _] in storage:
                                of.write(a + ' ')
                            of.write('\n')
                            for [_, a, _] in storage:
                                of.write(a + ' ')
                            of.write('\n')
                            for [_, _, a] in storage:
                                of.write(a + ' ')
                            of.write('\n')
                            storage = []

                    elif wait_faces > 0:
                        ms = l.split()
                        of.write(ms[0] + ' ' + ms[1] + ' ' + ms[2] + '\n')
                        of.write(ms[2] + ' ' + ms[3] + ' ' + ms[0] + '\n')
                        wait_faces -= 1

                l = f.readline()

            of.close()
        f.close()

#---------------------------------------------------------------------------------------------------

def step_0_add_data(filename, ofilename):
    """
    Add zero data.

    Arguments:
        filename -- input file name,
        ofilename -- output file name.
    """

    is_title_str = lambda s: 'TITLE' in s
    is_variables_str = lambda s: 'VARIABLES' in s
    is_zone_str = lambda s: 'ZONE' in s
    is_meta_str = lambda s: is_title_str(s) or is_variables_str(s) or is_zone_str(s)

    with open(filename, 'r') as f:
        with open(ofilename, 'w') as of:

            # Write head to out file.
            of.write('# EXPORT MODE: CHECK_POINT\n')
            of.write('TITLE = "FE Surface Data ASCII"\n')
            of.write('VARIABLES = "X", "Y", "Z", "T", "Hw", "Hi", "Node_HTC", "Node_Beta", "Node_TauX", "Node_TauY", "Node_TauZ"\n')

            wait_nodes, wait_faces = 0, 0
            l = f.readline()
            while l:

                if is_meta_str(l):
                    if (wait_nodes > 0) or (wait_faces > 0):
                        raise Exception('internal error')
                    if is_title_str(l):
                        print('ORIGINAL : ', l)
                    elif is_variables_str(l):
                        print('ORIGINAL : ', l)
                        variables = l[l.index('=') + 1 : ].split(',')
                        variables_count = len(variables)
                        print('variables = ', variables, ' count = ', variables_count)
                    elif is_zone_str(l):
                        print('ORIGINAL : ', l)
                        fields = l[l.index('NODES=') : ].split()
                        print('fields = ', fields)
                        nodes_count, faces_count = int(fields[0][6:-1]), int(fields[1][9:-1])
                        print('nodes_count = ', nodes_count, ' faces_count = ', faces_count)
                        wait_nodes, wait_faces = variables_count, faces_count
                        of.write('ZONE T="TRIANGLES", NODES=%d, ELEMENTS=%d, DATAPACKING=BLOCK, ZONETYPE=FETRIANGLE, VARLOCATION=([4-6]=CELLCENTERED)\n'
                                 % (nodes_count, faces_count))
                else:
                    if wait_nodes > 0:
                        of.write(l)
                        l = f.readline()
                        of.write(l)
                        l = f.readline()
                        of.write(l)
                        for i in range(3):
                            for _ in range(faces_count):
                                of.write('0.0 ')
                            of.write('\n')
                        for i in range(5):
                            for _ in range(nodes_count):
                                of.write('0.0 ')
                            of.write('\n')
                        wait_nodes = 0
                    elif wait_faces > 0:
                        of.write(l)
                        wait_faces -= 1

                l = f.readline()

            of.close()
        f.close()

#---------------------------------------------------------------------------------------------------

def correct_grid_data(g, filename):
    """
    Correct data.

    Arguments:
        g -- grid,
        filename -- name of file.
    """

    g.Nodes.sort()

    with open(filename, 'r') as f:

        l = f.readline()
        i = 0
        while l:

            s = l.split()
            x, y, z, beta, sf1, sf2, sf3, flux = float(s[0]), float(s[1]), float(s[2]), \
                                                 float(s[3]),                           \
                                                 float(s[4]), float(s[5]), float(s[6]), \
                                                 float(s[7])

            (nearest_node, dist) = g.NearestNodeInSortedNodes(x, y, z)
            if dist < 0.001:
                #print(i, ' | ', x, y, z, ' : ', str(nearest_node), ' : ', dist)
                nearest_node.Beta = beta
                nearest_node.Tau.X = sf1
                nearest_node.Tau.Y = sf2
                nearest_node.Tau.Z = sf3
                nearest_node.HTC = flux

            l = f.readline()
            i += 1

            if i % 1000 == 0:
                print('CHECK : %d' % i)

        f.close()

#---------------------------------------------------------------------------------------------------
# Test.
#---------------------------------------------------------------------------------------------------

if __name__ == '__main__':

    # Open data file and extract needed data.

    # X, Y, Z, Beta.
    #extract_data_from_dat('../data/air_inlet/drop3d_H9500_M06_T30.dat',
    #                      '../data/air_inlet/drop3d_H9500_M06_T30_out.dat',
    #                      [1, 1, 1, 0, 0, 0, 0, 1])

    # X, Y, Z, SF1, SF2, SF2, Classical Heat Flux.
    #extract_data_from_dat('../data/air_inlet/fensap_H9500_M06_T30.dat',
    #                      '../data/air_inlet/fensap_H9500_M06_T30_out.dat',
    #                      [1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0])

    # Merge files.
    # X, Y, Z, Beta, SF1, SF2, SF2, Classical Heat Flux.
    #merge_files('../data/air_inlet/drop3d_H9500_M06_T30_out.dat',
    #            '../data/air_inlet/fensap_H9500_M06_T30_out.dat',
    #            '../data/air_inlet/drop3d_fensap_out.dat')

    # Bring step_0 file to required format.
    # Only X, Y, Z, fields. Other are zeroes.
    #step_0_cut_extra_variables('../data/air_inlet/step_0.dat',
    #                           '../data/air_inlet/step_0_1.dat')

    # Change data packing type from point to block.
    #step_0_pack_point_to_block('../data/air_inlet/step_0_1.dat',
    #                           '../data/air_inlet/step_0_2.dat')

    # Add data to tecplot.
    #step_0_add_data('../data/air_inlet/step_0_2.dat',
    #                '../data/air_inlet/step_0_3.dat')

    # Load grid and data and save it back.
    g = gsu.Grid()
    g.ImportFromTecplot('../data/air_inlet/air_inlet_small.dat')
    correct_grid_data(g, '../data/air_inlet/drop3d_fensap_out.dat')
    g.ExportToTecplot('../data/air_inlet/air_inlet.dat')

#---------------------------------------------------------------------------------------------------
