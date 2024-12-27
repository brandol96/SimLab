import numpy as np
import os
import shutil


def write_waveplot(path,orbitals,KPTs,WP_grid,WP_Box_View,periodic):
    # its interesting here to create a file reader object holding 
    # the waveplot_input from "ugly_solutions" folder
    
    if periodic:
        WP = open('/home/rbrandolt/python-libs/SimLab/ugly_solutions/waveplot_in_periodic.hsd', 'r').readlines()
    else:
        WP = open('/home/rbrandolt/python-libs/SimLab/ugly_solutions/waveplot_in_molecule.hsd', 'r').readlines()
    
    with open(f'{path}waveplot_in.hsd','w+') as inp:
        for line in WP:
            if 'PlottedLevels' in line:
                inp.write(f'  PlottedLevels = {int(orbitals)}\n')
            elif 'PlottedKpoints' in line:
                inp.write(f'  PlottedKpoints = 1 {KPTs}\n')
            elif 'NrOfPoints' in line:
                inp.write(f'  NrOfPoints = {WP_grid} {WP_grid} {WP_grid}\n')
            elif 'RepeatBox' in line:
                inp.write(f'  RepeatBox = {{ {WP_Box_View[0]} {WP_Box_View[1]} {WP_Box_View[2]} }}\n')
            else:
                inp.write(line)

def run(Homo, Lumo, opt_path, orb_path, homo_max_kpt, lumo_min_kpt, WP_grid, WP_Box_View, periodic):
    i = 0
    N = len(Homo) - 1
    current_path = os.getcwd()

    for h in Homo:
        try:
            os.mkdir(f'{orb_path}homo-{N - i}')
        except:
            pass
        # HOMO
        print(f'current orbital: homo-{N - i}, {h}')
        shutil.copy(f'{opt_path}detailed.xml', f'{orb_path}homo-{N - i}{os.sep}detailed.xml')
        shutil.copy(f'{opt_path}eigenvec.bin', f'{orb_path}homo-{N - i}{os.sep}eigenvec.bin')
        write_waveplot(f'{orb_path}homo-{N - i}{os.sep}', h, homo_max_kpt, WP_grid, WP_Box_View, periodic)

        os.chdir(f'{orb_path}homo-{N - i}')
        os.system('waveplot > waveplot.out')
        os.chdir(current_path)
        i += 1

    i = 0
    N = len(Lumo) - 1
    for l in Lumo:
        try:
            os.mkdir(f'{orb_path}lumo+{i}')
        except:
            pass
        # LUMO
        print(f'current orbital: lumo+{i}, {l}')
        shutil.copy(f'{opt_path}detailed.xml', f'{orb_path}lumo+{i}{os.sep}detailed.xml')
        shutil.copy(f'{opt_path}eigenvec.bin', f'{orb_path}lumo+{i}{os.sep}eigenvec.bin')
        write_waveplot(f'{orb_path}lumo+{i}{os.sep}', l, lumo_min_kpt, WP_grid, WP_Box_View, periodic)

        os.chdir(f'{orb_path}lumo+{i}')
        os.system('waveplot > waveplot.out')
        os.chdir(current_path)
        i += 1
    print('')
