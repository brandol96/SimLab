import numpy as np
import os
import shutil


def write_waveplot(path,orbitals):
    with open('/home/rbrandolt/python-libs/SimLab/ugly_solutions/waveplot_in.hsd','r') as WP:
        with open(f'{path}waveplot_in.hsd','w+') as inp:
            for line in WP:
                if 'PlottedLevels' in line:
                    inp.write(f'  PlottedLevels = {int(orbitals)}\n')
                else:
                    inp.write(line)

def run(Homo, Lumo, opt_path, orb_path):
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
        write_waveplot(f'{orb_path}homo-{N - i}{os.sep}', h)

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
        write_waveplot(f'{orb_path}lumo+{i}{os.sep}', l)

        os.chdir(f'{orb_path}lumo+{i}')
        os.system('waveplot > waveplot.out')
        os.chdir(current_path)
        i += 1
    print('')