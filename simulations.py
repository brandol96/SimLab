# here are the set of stuff that a simulation needs done before starting
import shutil
import os
from ase.io import read
from ase.io import write


def start_sim(mol, mol_name, out_path, **kwargs):
    simulation = kwargs.get('simulation')
    method = kwargs.get('method')

    if method == 'DFTB':
        if simulation == 'optimize':
            from SimLab.calculator import fetch_dftb_calc
            print(f'{method} optimization for {mol_name}')
            pbc = mol.get_pbc()
            if True in pbc:
                print('Some direction has pbc, set "lattice_opt = True" to optimize the unit cell')
                dftb = fetch_dftb_calc(mol, cluster=False, **kwargs)
            else:
                print('No direction has pbc, DFTB will NOT perform lattice optimization')
                dftb = fetch_dftb_calc(mol, cluster=True, **kwargs)
            # run optimization through DFTB+ implemented routines
            mol.set_calculator(dftb)
            dftb.calculate(mol)

            mol = read(f'geo_end.gen')
            mol.center()
            write(f'{method}_{mol_name}_end.traj', mol)
        elif simulation == 'optical':
            from SimLab.analysis import optical
            print(f'{method} optimization for {mol_name}')
            curr_ase_dftb_command = os.environ["ASE_DFTB_COMMAND"]
            os.environ["ASE_DFTB_COMMAND"] = "dftb+ | tee PREFIX.out"

            max_SCC = kwargs.get('max_SCC')
            max_SCC_steps = kwargs.get('max_SCC_steps')
            fermi_filling = kwargs.get('fermi_filling')
            total_time = kwargs.get('total_time')
            time_step = kwargs.get('time_step')
            field_strength = kwargs.get('field_strength')
            directions = kwargs.get('directions')
            laser = kwargs.get('laser')
            laser_energy = kwargs.get('laser_energy')
            n_points = kwargs.get('n_points')

            pbc = mol.get_pbc()
            if True in pbc:
                print('Some direction has pbc, the molecule is NOT valid! ( Yet :o ) \n\n')
            else:
                print('No direction has pbc, the molecule is valid!')
                opt_out_path = f'optimize_{method}_{mol_name}' + os.sep
                shutil.copyfile(f'{opt_out_path}charges.bin', f'{os.getcwd()}{os.sep}charges.bin')

                if not laser:
                    # calculation using kick expects direction to be a string containing
                    # the desired directions 'XYZ', 'XY', 'ZY', etc...
                    for direction in directions:
                        optical.run_kick(mol, max_SCC, max_SCC_steps, fermi_filling,
                                         total_time, time_step, field_strength, n_points, direction)
                else:
                    # calculation using kick expects direction a list of direction intensities
                    # such as [0.0, 1.0 ,0.0] -> laser only in Y direction [0.5 0.0 0.5] -> half X and Z directions
                    optical.run_laser(mol, max_SCC, max_SCC_steps, fermi_filling,
                                      total_time, time_step, laser_energy, field_strength, n_points, directions)
                print('\n\n')
            os.environ["ASE_DFTB_COMMAND"] = curr_ase_dftb_command
        else:
            print(f'{simulation} is not possible with {method}')


def start_view(mol, mol_name, out_path, **kwargs):
    simulation = kwargs.get('simulation')
    method = kwargs.get('method')
    if method == 'DFTB':
        if simulation == 'optical':
            from SimLab.view import optical
            print(f'{method} optimization for {mol_name}')
            curr_ase_dftb_command = os.environ["ASE_DFTB_COMMAND"]
            os.environ["ASE_DFTB_COMMAND"] = "dftb+ | tee PREFIX.out"
            field_strength = kwargs.get('field_strength')
            method = kwargs.get('method')
            directions = kwargs.get('directions')
            laser = kwargs.get('laser')
            interactive_plot = kwargs.get('interactive_plot')
            fourier_damp = kwargs.get('fourier_damp')
            laser_energy = kwargs.get('laser_energy')

            pbc = mol.get_pbc()
            if True in pbc:
                print('Some direction has pbc, the molecule is NOT valid! ( Yet :o ) \n\n')
            else:
                print('No direction has pbc, the molecule is valid!')

                if not laser:
                    # the output info from kick optical simulation is the full spectra
                    # as well as the location of possible absorption peaks
                    optical.run_kick(method, out_path, mol_name,
                                     interactive_plot, directions, laser, fourier_damp, field_strength)
                else:
                    # the output info from laser optical simulation is the polarization response
                    # of the material to the laser excitation and
                    # TODO: decide if charge distribution or potential will be the relevant output
                    optical.run_laser(method, out_path, mol_name,
                                      interactive_plot, directions, laser, fourier_damp, field_strength)
                print('\n\n')
            os.environ["ASE_DFTB_COMMAND"] = curr_ase_dftb_command
