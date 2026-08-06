import numpy as np
import os
import shutil
def run_waveplot_parallelism(OMP_threads,MPI_cores,verbosity):
    if MPI_cores != 1:
        dftb_mpi_bin = os.environ["DFTB_MPI_BIN"]
        dftb_mpi_lib = os.environ["DFTB_MPI_LIB"]
        inline_env = f"LD_LIBRARY_PATH={dftb_mpi_lib}"
        if verbosity > 2:
            os.system(f'{inline_env} {dftb_mpi_bin}waveplot | tee waveplot.out')
        else:
            os.system(f'{inline_env} {dftb_mpi_bin}waveplot > waveplot.out')
        return
    else:
        os.environ["OMP_PROC_BIND"] = 'true'
        os.environ["OMP_PLACES"] = 'cores'
        os.environ["OPENBLAS_NUM_THREADS"] = str(OMP_threads)
        os.environ["OMP_NUM_THREADS"] = str(OMP_threads)
        dftb_omp_bin = os.environ["DFTB_OMP_BIN"]
        dftb_omp_lib = os.environ["DFTB_OMP_LIB"]
        inline_env = f"LD_LIBRARY_PATH={dftb_omp_lib}"
        if verbosity > 2:
            os.system(f'{inline_env} {dftb_omp_bin}waveplot | tee waveplot.out')
        else:
            os.system(f'{inline_env} {dftb_omp_bin}waveplot > waveplot.out')
        return


def write_waveplot(path,orbitals,KPTs,WP_grid,WP_Box_View,periodic,SKfiles):
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
            elif '<<+' in line:
                inp.write(f'<<+ "{SKfiles}wfc.hsd"\n')
            else:
                inp.write(line)

def run(Homo, Lumo, opt_path, orb_path, homo_max_kpt, lumo_min_kpt,
        WP_grid, WP_Box_View, periodic, SKfiles,
        OMP_threads,MPI_cores,verbosity):
    i = 0
    N = len(Homo) - 1
    current_path = os.getcwd()

    #TODO: there is a possibility that HOMO-n <= 0. This orbital does not make sense and should be tested for

    for h in Homo:
        try:
            os.mkdir(f'{orb_path}homo-{N - i}')
        except:
            pass
        # HOMO
        print(f'current orbital: homo-{N - i}, {h}')
        shutil.copy(f'{opt_path}detailed.xml', f'{orb_path}homo-{N - i}{os.sep}detailed.xml')
        shutil.copy(f'{opt_path}eigenvec.bin', f'{orb_path}homo-{N - i}{os.sep}eigenvec.bin')
        write_waveplot(f'{orb_path}homo-{N - i}{os.sep}', h, homo_max_kpt, WP_grid, WP_Box_View, periodic,SKfiles)

        os.chdir(f'{orb_path}homo-{N - i}')
        run_waveplot_parallelism(OMP_threads,MPI_cores,verbosity)
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
        write_waveplot(f'{orb_path}lumo+{i}{os.sep}', l, lumo_min_kpt, WP_grid, WP_Box_View, periodic,SKfiles)

        os.chdir(f'{orb_path}lumo+{i}')
        run_waveplot_parallelism(OMP_threads,MPI_cores,verbosity)
        os.chdir(current_path)
        i += 1
    print('')
