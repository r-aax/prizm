# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 15:45:02 2019

@author: Rybakov
"""

from functools import reduce
import lst

#---------------------------------------------------------------------------------------------------
# Functions.
#---------------------------------------------------------------------------------------------------

def is_title_str(s):
    """
    Check if string is title.

    Arguments:
        s -- string.

    Result:
        True -- is string is title,
        False -- otherwise.
    """

    return 'TITLE' in s

#---------------------------------------------------------------------------------------------------

def is_variables_str(s):
    """
    Check if string is variables.

    Arguments:
        s -- string.

    Result:
        True -- is string is variables,
        False -- otherwise.
    """

    return 'VARIABLES' in s

#---------------------------------------------------------------------------------------------------

def is_zone_str(s):
    """
    Check if string is zone.

    Arguments:
        s -- string.

    Result:
        True -- is string is zone,
        False -- otherwise.
    """

    return 'ZONE' in s

#---------------------------------------------------------------------------------------------------

def is_meta_str(s):
    """
    Check if string is meta.

    Arguments:
        s -- string.

    Result:
        True -- is string is meta,
        False -- otherwise.
    """

    return is_title_str(s) or is_variables_str(s) or is_zone_str(s)

#---------------------------------------------------------------------------------------------------

def unpack_variables(s):
    """
    Extract variables from string.

    Arguments:
        s -- string.

    Result:
        Variables names list.
    """

    return s[s.index('=') + 1 : ].split(',')

#---------------------------------------------------------------------------------------------------

def pack_variables(vs):
    """
    Pack variables to string.

    Arguments:
        vs -- variables names.

    Result:
        Variables string.
    """

    vs = [''.join(list(filter(lambda x: not x.isspace(), v))) for v in vs]

    return 'VARIABLES=' + reduce(lambda x, y: x + ', ' + y, vs)

#---------------------------------------------------------------------------------------------------

def pack_str_array(a):
    """
    Pack strings array.

    Arguments:
        a -- array.

    Result:
        String.
    """

    return reduce(lambda x, y: x + ' ' + y, a)

#---------------------------------------------------------------------------------------------------

def get_zonetype_nodes_elements(s):
    """
    Get zonetype, nodes count and elements count from string.

    Arguments:
        s -- string.

    Result:
        Tuple (zonetype, nodes count, elements count).
    """

    fields = s.split('=')
    fields = fields[1] + fields[4] + fields[5]
    fields = fields.split()
    field_z = fields[0][:-1]
    field_n = fields[2]
    field_e = fields[5]

    return (field_z, int(field_n), int(field_e))

#---------------------------------------------------------------------------------------------------

def select_variables(filename, ofilename, mask):
    """
    Select variables from *.dat file.

    Arguments:
        filename -- input file,
        ofilename -- output file,
        mask -- variables mask.
    """

    print('select_variables : ', filename, ' -> ', ofilename)

    with open(filename, 'r') as f:
        with open(ofilename, 'w') as of:

            # Loop for all input file lines.
            zonetype= ''
            wait_nodes, wait_elements = 0, 0
            l = f.readline()
            while l:

                if is_meta_str(l):
                    if (wait_nodes > 0) or (wait_elements > 0):
                        raise Exception('internal error')
                    if is_title_str(l):
                        print('ORIGINAL TITLE : ', l)
                        of.write(l)
                    elif is_variables_str(l):
                        print('ORIGINAL VARIABLES : ', l)
                        ol = pack_variables(lst.mask(unpack_variables(l), mask))
                        print('NEW VARIABLES : ', ol)
                        of.write(ol + '\n')
                    elif is_zone_str(l):
                        print('ORIGINAL ZONE : ', l)
                        (zonetype, nodes_count, elements_count) = get_zonetype_nodes_elements(l)
                        print('ZONE STAT : ', zonetype, nodes_count, elements_count)
                        if zonetype != 'FEBRICK':
                            of.write(l)
                            wait_nodes, wait_elements = nodes_count, elements_count
                else:
                    if wait_nodes > 0:
                        ms = lst.mask(l.split(), mask)
                        if ms != []:
                            out = reduce(lambda a, b: a + ' ' + b, ms) + '\n'
                            of.write(out)
                        wait_nodes -= 1
                    elif wait_elements > 0:
                        of.write(l)
                        wait_elements -= 1

                l = f.readline()

    f.close()
    of.close()

#---------------------------------------------------------------------------------------------------

def filter_surfaces(filename, ofilename, indices):
    """
    Filter surfaces.

    Arguments:
        filename -- input file,
        ofilename -- output file,
        indices -- indices list
    """

    print('filter_surfaces : ', filename, ' -> ', ofilename)

    with open(filename, 'r') as f:
        with open(ofilename, 'w') as of:

            # Loop for all input file lines.
            zone_index = 0
            wait_nodes, wait_elements = 0, 0
            l = f.readline()
            while l:

                if is_meta_str(l):
                    if (wait_nodes > 0) or (wait_elements > 0):
                        raise Exception('internal error')
                    if is_title_str(l):
                        of.write(l)
                    elif is_variables_str(l):
                        of.write(l)
                    elif is_zone_str(l):
                        #print('ORIGINAL ZONE : ', l)
                        (zonetype, nodes_count, elements_count) = get_zonetype_nodes_elements(l)
                        #print('ZONE STAT : ', zonetype, nodes_count, elements_count)
                        if (indices == None) or (zone_index in indices):
                            print('INDEX : ', zone_index, '+')
                            of.write(l)
                            wait_nodes, wait_elements = nodes_count, elements_count
                        else:
                            print('INDEX : ', zone_index, '-')
                        zone_index += 1
                else:
                    if wait_nodes > 0:
                        of.write(l)
                        wait_nodes -= 1
                    elif wait_elements > 0:
                        of.write(l)
                        wait_elements -= 1

                l = f.readline()

    f.close()
    of.close()

#---------------------------------------------------------------------------------------------------

def merge_drop_fens(drop_filename, fens_filename, ofilename):
    """
    Merge drop and fens into compressor.dat.

    Arguments:
        drop_filename -- drop filename,
        fens_filename -- fens filename,
        ofilename -- output filename.
    """

    print('merge_drop_fens : ', drop_filename, ', ', fens_filename, ' -> ', ofilename)

    with open(drop_filename, 'r') as df:
        with open(fens_filename, 'r') as ff:
            with open(ofilename, 'w') as of:

                # Loop for all input file lines.
                wait_nodes, wait_elements = 0, 0
                dl, fl = df.readline(), ff.readline()
                while dl:

                    if is_meta_str(dl):
                        if (wait_nodes > 0) or (wait_elements > 0):
                            raise Exception('internal error')
                        if is_title_str(dl):
                            of.write(dl)
                        elif is_variables_str(dl):
                            of.write('VARIABLES = "X", "Y", "Z", ' \
                                     + '"Node_Flux", "Node_Beta", ' \
                                     + '"Node_TauX", "Node_TauY", "Node_TauZ"' + '\n')
                        elif is_zone_str(dl):
                            (zonetype, nodes_count, elements_count) = get_zonetype_nodes_elements(dl)
                            of.write(dl)
                            wait_nodes, wait_elements = nodes_count, elements_count
                    else:
                        if wait_nodes > 0:
                            dfields = dl.split()
                            ffields = fl.split()
                            result_line = dfields[0] + ' ' + dfields[1] + ' ' + dfields[2] + ' ' \
                                          + ffields[6] + ' ' + dfields[3] + ' ' \
                                          + ffields[3] + ' ' + ffields[4] + ' ' + ffields[5] + '\n'
                            of.write(result_line)
                            wait_nodes -= 1
                        elif wait_elements > 0:
                            of.write(dl)
                            wait_elements -= 1

                    dl, fl = df.readline(), ff.readline()


    df.close()
    ff.close()
    of.close()

#---------------------------------------------------------------------------------------------------

def calculate_htc(filename, ofilename):
    """
    Calculate HTC from Flux.

    Arguments:
        filename -- name of input file,
        ofilename -- name of output file.
    """

    print('calculate_htc : ', filename, ' -> ', ofilename)

    with open(filename, 'r') as f:
        with open(ofilename, 'w') as of:

            # Loop for all input file lines.
            wait_nodes, wait_elements = 0, 0
            l = f.readline()
            while l:

                if is_meta_str(l):
                    if (wait_nodes > 0) or (wait_elements > 0):
                        raise Exception('internal error')
                    if is_title_str(l):
                        of.write(l)
                    elif is_variables_str(l):
                        of.write('VARIABLES = "X", "Y", "Z", ' \
                                 + '"Node_HTC", "Node_Beta", ' \
                                 + '"Node_TauX", "Node_TauY", "Node_TauZ"' + '\n')
                    elif is_zone_str(l):
                        (zonetype, nodes_count, elements_count) = get_zonetype_nodes_elements(l)
                        wait_nodes, wait_elements = nodes_count, elements_count
                        of.write(l)
                else:
                    if wait_nodes > 0:
                        fields = l.split()
                        flux = float(fields[3])

                        # Calc HTC.
                        temperature = 255.15 - 273.15
                        velocity = 184.0
                        cp = 1009.0
                        t_adiabatic = temperature \
                                      + 0.5 * velocity * velocity / cp
                        t_inc = 10.0
                        t_wall = t_adiabatic + t_inc
                        dtemp = t_wall - temperature
                        htc = abs(flux / dtemp)

                        res = fields[0] + ' ' + fields[1] + ' ' + fields[2] + ' ' \
                              + str(htc) + ' ' + fields[4] + ' ' \
                              + fields[5] + ' ' + fields[6] + ' ' + fields[7] + '\n'
                        of.write(res)
                        wait_nodes -= 1
                    elif wait_elements > 0:
                        of.write(l)
                        wait_elements -= 1

                l = f.readline()

    f.close()
    of.close()

#---------------------------------------------------------------------------------------------------

def convert_to_block_pack(filename, ofilename):
    """
    Convert grid to block packing mode.

    # EXPORT MODE: CHECK_POINT
    TITLE="FE Surface Data ASCII"
    VARIABLES="X", "Y", "Z", "T", "Hw", "Hi", "Node_HTC", "Node_Beta", "Node_TauX", "Node_TauY", "Node_TauZ"
    ZONE T="TRIANGLES"
    NODES=396
    ELEMENTS=394
    DATAPACKING=BLOCK
    ZONETYPE=FETRIANGLE
    VARLOCATION=([4-6]=CELLCENTERED)

    Arguments:
        filename -- name of file,
        ofilename -- name of output file.
    """

    print('convert_to_block_pack : ', filename, ' -> ', ofilename)

    with open(filename, 'r') as f:
        with open(ofilename, 'w') as of:

            # Loop for all input file lines.
            nodes_count, elements_count = 0, 0
            wait_nodes, wait_elements = 0, 0
            l = f.readline()
            while l:

                if is_meta_str(l):
                    if (wait_nodes > 0) or (wait_elements > 0):
                        raise Exception('internal error')
                    if is_title_str(l):
                        of.write('# EXPORT MODE: CHECK_POINT\n')
                        of.write('TITLE="FE Surface Data ASCII"\n')
                    elif is_variables_str(l):
                        of.write('VARIABLES="X", "Y", "Z", ' \
                                 + '"T", "Hw", "Hi", ' \
                                 + '"Node_HTC", "Node_Beta", ' \
                                 + '"Node_TauX", "Node_TauY", "Node_TauZ"' + '\n')
                    elif is_zone_str(l):
                        print('ZONE')
                        (zonetype, nodes_count, elements_count) = get_zonetype_nodes_elements(l)
                        of.write('ZONE T="TRIANGLES"\n')
                        of.write('NODES=%d\n' % nodes_count)
                        of.write('ELEMENTS=%d\n' % (2 * elements_count))
                        of.write('DATAPACKING=BLOCK\n')
                        of.write('ZONETYPE=FETRIANGLE\n')
                        of.write('VARLOCATION=([4-6]=CELLCENTERED)\n')
                        wait_nodes, wait_elements = nodes_count, elements_count
                        xs = [''] * nodes_count
                        ys = [''] * nodes_count
                        zs = [''] * nodes_count
                        htcs = [''] * nodes_count
                        betas = [''] * nodes_count
                        tauxs = [''] * nodes_count
                        tauys = [''] * nodes_count
                        tauzs = [''] * nodes_count
                else:
                    if wait_nodes > 0:
                        fields = l.split()
                        xs[nodes_count - wait_nodes] = fields[0]
                        ys[nodes_count - wait_nodes] = fields[1]
                        zs[nodes_count - wait_nodes] = fields[2]
                        htcs[nodes_count - wait_nodes] = fields[3]
                        betas[nodes_count - wait_nodes] = fields[4]
                        tauxs[nodes_count - wait_nodes] = fields[5]
                        tauys[nodes_count - wait_nodes] = fields[6]
                        tauzs[nodes_count - wait_nodes] = fields[7]
                        wait_nodes -= 1

                        if wait_nodes == 0:
                            
                            # Final print.
                            of.write(pack_str_array(xs) + '\n')
                            of.write(pack_str_array(ys) + '\n')
                            of.write(pack_str_array(zs) + '\n')
                            of.write(pack_str_array(['0.0'] * (2 * elements_count)) + '\n')
                            of.write(pack_str_array(['0.0'] * (2 * elements_count)) + '\n')
                            of.write(pack_str_array(['0.0'] * (2 * elements_count)) + '\n')
                            of.write(pack_str_array(htcs) + '\n')
                            of.write(pack_str_array(betas) + '\n')
                            of.write(pack_str_array(tauxs) + '\n')
                            of.write(pack_str_array(tauys) + '\n')
                            of.write(pack_str_array(tauzs) + '\n')

                    elif wait_elements > 0:
                        fields = l.split()
                        of.write(fields[0] + ' ' + fields[2] + ' ' + fields[1] + '\n')
                        of.write(fields[0] + ' ' + fields[3] + ' ' + fields[2] + '\n')
                        wait_elements -= 1

                l = f.readline()

    f.close()
    of.close()

#---------------------------------------------------------------------------------------------------
# Test.
#---------------------------------------------------------------------------------------------------

if __name__ == '__main__':

    print('prepare_cylinder_for_crystal')

    stage = 0

    # Этап 1.
    # Из файла droplet выбираем только нужные нам переменные.
    # Сохраняем структуру multizone.
    # Выбираем только "X", "Y", "Z",
    #                 4 переменные пропускаем,
    #                 "Collection efficiency-Droplet".
    if stage == 1:
        filename = '../data/cylinder/droplet.dat'
        ofilename = '../data/cylinder/droplet_surfaces.dat'
        select_variables(filename, ofilename, [1, 1, 1, 0, 0, 0, 0, 1])

    # Этап 2.
    # Из файла fensap выбираем только нужные нам переменные.
    # Сохраняем структуру multizone.
    # Выбираем только "X", "Y", "Z",
    #                 12 переменных пропускаем,
    #                 3 переменные shear-stress,
    #                 "Classical heat flux",
    #                 1 переменную пропускаем.
    if stage == 2:
        filename = '../data/cylinder/fensap.dat'
        ofilename = '../data/cylinder/fensap_surfaces.dat'
        select_variables(filename, ofilename,
                         [1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0])

    # Этап 3.
    # Фильтрация поверхностей droplet_surfaces.
    # 1.
    if stage == 3:
        filter_surfaces('../data/cylinder/droplet_surfaces.dat',
                        '../data/cylinder/droplet_filtered_surfaces.dat',
                        [1])

    # Этап 4.
    # Фильтрация поверхностей fensap_surfaces.
    # 1.
    if stage == 4:
        filter_surfaces('../data/cylinder/fensap_surfaces.dat',
                        '../data/cylinder/fensap_filtered_surfaces.dat',
                        [1])

    # Этап 5.
    # Объединение зон в файлы drop и fens.
    # Выпоняется в bash.
    pass

    # Этап 6.
    # Слияние drop и fens с правильным порядком переменных.
    if stage == 6:
        merge_drop_fens('../data/cylinder/droplet_filtered_surfaces.dat',
                        '../data/cylinder/fensap_filtered_surfaces.dat',
                        '../data/cylinder/cylinder_point_pack.dat')

    # Этап 7.
    # Пересчет HTC из потока тепла.
    if stage == 7:
        calculate_htc('../data/cylinder/cylinder_point_pack.dat',
                      '../data/cylinder/cylinder_htc.dat')

    # Этап 8.
    # Переработка в формат crystal.
    if stage == 8:
        convert_to_block_pack('../data/cylinder/cylinder_htc.dat',
                              '../data/cylinder/cylinder.dat')
