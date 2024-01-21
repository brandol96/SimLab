# here are the set of stuff that a simulation needs done before starting
import shutil
import os
from ase.io import read
from ase.io import write


# to choose a simulation the varaible 'simulation' msut be fileld along with any
# required variables, any value left undeclared will be assumed to be standard
# notice not all standard values work with all possible simulations
# available so far:
# -> optimize
# -> molecular_dynamics
# -> bands
# -> optical
# -----> kick
# -----> laser
# -----> casida

def start_sim(mol, mol_name, out_path, **kwargs):
    simulation = kwargs.get('simulation')
    method = kwargs.get('method')

    if method == 'DFTB':
        SKFiles = kwargs.get('SKFiles')
        os.environ["DFTB_PREFIX"] = f"/home/rbrandolt/dftb-23.1/skfiles/{SKFiles}"
        if simulation == 'optimize':
            from SimLab.calculator import fetch_dftb_calc
            print(f'Run: {method} optimization for {mol_name}')
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

            try:
                mol = read(f'geo_end.gen')
                mol.center()
                write(f'{method}_{mol_name}_end.traj', mol)
            except:
                print('First step did not converge SCC!')

        elif simulation == 'molecular_dynamics':
            from SimLab.analysis import molecular_dynamics

            kpts = kwargs.get('kpts')
            thermostat = kwargs.get('thermostat')
            temp_profile = kwargs.get('temp_profile')
            time_step = kwargs.get('time_step')
            SCC = kwargs.get('SCC')
            max_SCC = kwargs.get('max_SCC')
            max_SCC_steps = kwargs.get('max_SCC_steps')
            fermi_filling = kwargs.get('fermi_filling')

            print(f'Run: {method} molecular dynamics for {mol_name}')
            molecular_dynamics.run(mol, mol_name, kpts, thermostat, temp_profile, time_step,
                                   SCC, max_SCC, max_SCC_steps, fermi_filling)

        elif simulation == 'bands':
            print(f'{method} band structure for {mol_name}')
            opt_out_path = f'optimize_{method}_{mol_name}{os.sep}'
            shutil.copyfile(f'{opt_out_path}charges.bin', f'{os.getcwd()}{os.sep}charges.bin')
            pbc = mol.get_pbc()

            if True in pbc:
                from SimLab.utils import path_dftb
                from SimLab.analysis import bands
                BZ_path = kwargs.get('BZ_path')
                BZ_step = kwargs.get('BZ_step')
                verbosity = kwargs.get('verbosity')
                sampling = path_dftb(BZ_path, BZ_step, mol, verbosity, False)
                print(f'\nSome direction has pbc, I\'ll use the provided path: {BZ_path} and step: {BZ_step}')
                print(f'Path generated:\n\n{sampling}')
                print(f'Use verbosity >= 3 for more details...')
                bands.run(mol, mol_name, sampling)

            print('\n\n')

        elif simulation == 'optical':
            from SimLab.analysis import optical
            max_SCC = kwargs.get('max_SCC')
            max_SCC_steps = kwargs.get('max_SCC_steps')
            fermi_filling = kwargs.get('fermi_filling')
            total_time = kwargs.get('total_time')
            time_step = kwargs.get('time_step')
            field_strength = kwargs.get('field_strength')
            directions = kwargs.get('directions')
            sim_type = kwargs.get('sim_type')
            laser_energy = kwargs.get('laser_energy')
            n_points = kwargs.get('n_points')
            n_excitations = kwargs.get('n_excitations')
            cutoff_energy = kwargs.get('cutoff_energy')
            cutoff_oscillator = kwargs.get('cutoff_osc')

            print(f'Run: {method} {sim_type} optical absorption for {mol_name}')
            pbc = mol.get_pbc()
            if True in pbc:
                print('Some direction has pbc, the molecule is NOT valid! ( Yet :o ) \n\n')
            else:
                print('No direction has pbc, the molecule is valid!')
                opt_out_path = f'optimize_{method}_{mol_name}' + os.sep
                shutil.copyfile(f'{opt_out_path}charges.bin', f'{os.getcwd()}{os.sep}charges.bin')
                if sim_type == 'kick':
                    # calculation using kick expects direction to be a string containing
                    # the desired directions 'XYZ', 'XY', 'ZY', etc...
                    for direction in directions:
                        print(f'current direction: {direction}')
                        optical.run_kick(mol, max_SCC, max_SCC_steps, fermi_filling,
                                         total_time, time_step, field_strength, n_points, direction)
                elif sim_type == 'laser':
                    print(f'current direction: {directions}')
                    # calculation using kick expects direction a list of direction intensities
                    # such as [0.0, 1.0 ,0.0] -> laser only in Y direction [0.5 0.0 0.5] -> half X and Z directions
                    optical.run_laser(mol, max_SCC, max_SCC_steps, fermi_filling,
                                      total_time, time_step, laser_energy, field_strength, n_points, directions)
                elif sim_type == 'casida':
                    print(f'Running Casida stuff!')
                    optical.run_casida(mol, max_SCC, max_SCC_steps, fermi_filling,
                                       n_excitations, cutoff_energy, cutoff_oscillator)

                print('\n\n')
        else:
            print(f'Run: {simulation} is not possible with {method}')


def start_view(mol, mol_name, out_path, **kwargs):
    simulation = kwargs.get('simulation')
    method = kwargs.get('method')
    interactive_plot = kwargs.get('interactive_plot')
    if method == 'DFTB':
        if simulation == 'optimize':
            from SimLab.view import DOS
            pbc = mol.get_pbc()
            dos_range = kwargs.get('dos_range')
            if True in pbc:
                print('Some direction has pbc, set "simulation = bands" to view band structure and DOS')
            else:
                print('No direction has pbc, simple DOS')
                DOS.run_dftb(mol, mol_name, out_path, interactive_plot, dos_range)

        if simulation == 'bands':
            from SimLab.view import bands
            print(f'View: {method} band structure for {mol_name}')
            pbc = mol.get_pbc()
            if True in pbc:
                BZ_path = kwargs.get('BZ_path')
                BZ_step = kwargs.get('BZ_step')
                bands_zoom = kwargs.get('bands_zoom')
                verbosity = kwargs.get('verbosity')
                print(f'\nSome direction has pbc, I\'ll use the provided path: {BZ_path} and step: {BZ_step}')
                bands.run_dftb(mol, mol_name, out_path, bands_zoom, BZ_path, BZ_step,
                               interactive_plot,verbosity)
            else:
                print('No direction has pbc, the molecule is invalid you silly!')

        if simulation == 'optical':
            from SimLab.view import optical
            sim_type = kwargs.get('sim_type')
            field_strength = kwargs.get('field_strength')
            method = kwargs.get('method')
            directions = kwargs.get('directions')
            interactive_plot = kwargs.get('interactive_plot')
            fourier_damp = kwargs.get('fourier_damp')
            cutoff_OscStr = kwargs.get('cutoff_OscStr')
            cutoff_weight = kwargs.get('cutoff_weight')
            energy_upper_plot = kwargs.get('energy_upper_plot')
            print(f'{method} {sim_type} optical absorption for {mol_name}')

            pbc = mol.get_pbc()
            if True in pbc:
                print('Some direction has pbc, the molecule is NOT valid! ( Yet :o ) \n\n')
            else:
                print('No direction has pbc, the molecule is valid!')

                if sim_type == 'kick':
                    # the output info from kick optical simulation is the full spectra
                    # as well as the location of possible absorption peaks
                    optical.run_kick(method, out_path, mol_name,
                                     interactive_plot, directions, fourier_damp, field_strength)
                elif sim_type == 'laser':
                    # the output info from laser optical simulation is the polarization response
                    # of the material to the laser excitation and
                    # TODO: decide if charge distribution or potential will be the relevant output
                    optical.run_laser(method, out_path, mol_name, interactive_plot)
                elif sim_type == 'casida':
                    # TODO: remove 'kick' requirement
                    print('Casida')
                    # this method depends on 'kick' optical to be already placed
                    # I know this is not really  required, but will do for now
                    #optical.run_casida(method, out_path, mol_name, cutoff_OscStr, cutoff_weight,
                    #                   interactive_plot, energy_upper_plot)
                    optical.run_casida_pure(method, out_path, mol_name,
                                            cutoff_OscStr, cutoff_weight, energy_upper_plot,
                                            interactive_plot)
                print('\n\n')

        elif simulation == 'orbitals':
            import numpy as np
            from SimLab.view import orbitals
            from SimLab.utils import read_fermi_levels_dftb

            # get required variables
            opt_out_path = f'optimize_{method}_{mol_name}' + os.sep
            target_orbirals = kwargs.get('target_orbirals')

            homo, lumo, gap, fermi_e = read_fermi_levels_dftb(opt_out_path, mol_name, verbosity=0)

            # setup target folder for orbitals
            H_0 = homo[1]
            L_0 = lumo[1]

            i = target_orbirals
            homo_list = np.arange(H_0 - i, H_0 + 1, 1)
            lumo_list = np.arange(L_0, L_0 + i + 1, 1)
            orb_path = f'orbitals_DFTB_{mol_name}{os.sep}'
            try:
                os.mkdir(orb_path)
            except:
                pass

            # start sim
            print(f'Run: {method} calculation for {mol_name}')
            orbitals.run(homo_list, lumo_list, opt_out_path, orb_path)