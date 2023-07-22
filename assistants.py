# This contains the basic class of Cody, used to control all of the performed simulations
# every folder file juggling should be places here to be done and written once.
import os


# noinspection PyUnresolvedReferences
class Cody:
    # This method sets required initial variables
    # I'll put the simulation parameters here
    def __init__(self, **kwargs):
        # internal setup parameters
        self.voice = kwargs.get('voice', False)
        print("\n\nHello, I am Cody!\n\n")
        if self.voice:
            os.system('spd-say "Hello, I am Cody!"')

        threads = kwargs.get('threads', 1)
        os.environ["OMP_NUM_THREADS"] = str(threads)
        self.verbosity = kwargs.get('verbosity', 0)

        # simulation parameters
        self.parameters = dict(
            # generic stuff
            simulation=kwargs.get('simulation', 'Optimize'),
            view_only=kwargs.get('view_only', False),
            interactive_plot=kwargs.get('interactive_plot', False),
            method=kwargs.get('method', 'DFTB'),
            label=kwargs.get('label', 'dftb_output'),
            kpts=kwargs.get('kpts', (4, 4, 4)),
            lattice_opt=kwargs.get('lattice_opt', False),
            fix_angles=kwargs.get('fix_angles', False),
            fix_lengths=kwargs.get('fix_lengths', [False, False, False]),
            max_force=kwargs.get('max_force', 1E-4),
            time_step=kwargs.get('time_step', 0.005),

            # DFTB stuff
            max_driver_steps=kwargs.get('max_driver_steps', 10000),
            SCC=kwargs.get('SCC', True),
            max_SCC=kwargs.get('max_SCC', 1E-2),
            max_SCC_steps=kwargs.get('max_SCC_steps', 1000),
            fermi_filling=kwargs.get('fermi_filling', 0.0),
            use_LennardJones=kwargs.get('use_LennardJones', False),
            SKFiles=kwargs.get('SKFiles', 'Please supply a choice!'),

            # MD stuff
            thermostat=kwargs.get('thermostat', 'NVE'),
            temp_profile=kwargs.get('temp_profile', None),

            # band structure stuff
            BZ_path=kwargs.get('BZ_path', 'Please Supply a Path'),
            BZ_step=kwargs.get('BZ_step', 1E-2),
            bands_zoom=kwargs.get('bands_zoom', -1),

            # elastic constants stuff
            maxCauchyStrain=kwargs.get('maxCauchyStrain', 0.01),
            totalCauchySteps=kwargs.get('totalCauchySteps', 10),

            # optical absorption stuff
            sim_type=kwargs.get('sim_type', None), # I might use this parameter as generic
                                                   # to point at various methods of the same
                                                   # "kind" of simulation
            total_time=kwargs.get('total_time', 100),
            n_excitations=kwargs.get('n_excitations', 20),
            field_strength=kwargs.get('field_strength', 1E-3),
            directions=kwargs.get('directions', 'XYZ'),
            fourier_damp=kwargs.get('fourierDamp', 10),
            laser_energy=kwargs.get('laser_energy', 'PLEASE SUPPLY SOMETHING!'),
            cutoff_energy=kwargs.get('cutoff_energy', 30),
            cutoff_osc=kwargs.get('cutoff_osc', 0.000001),
            cutoff_OscStr=kwargs.get('cutoff_OscStr', 0.001),
            cutoff_weight=kwargs.get('cutoff_weight', 0.8),
            energy_upper_plot=kwargs.get('energy_upper_plot', 2.0),

            # coulomb potential
            n_points=kwargs.get('n_points', 1),
        )

        # TODO: check-up inconsistent parameters

    def change_instruction(self, par, new_value):
        self.parameters[par] = new_value
        # run re-check to ensure the user is not silly

    @staticmethod
    def fetch_molecule_list():
        input_list = os.listdir()
        mol_list = []
        for inputFile in input_list:
            if '.traj' in inputFile:
                mol_list.append(inputFile)
        return mol_list

    @staticmethod
    def clean_files(out_path):
        print('\n##### CLEANUP START ######\n')
        warning = True
        output_list = os.listdir()
        for outFile in output_list:
            keep = ('.traj' in outFile and '_end.traj' not in outFile) or \
                   '.py' in outFile or \
                   os.path.isdir(outFile) or \
                   'FermiLevels.out' == outFile or \
                   'effMass.out' == outFile

            if os.path.isdir(out_path):
                if warning:
                    print(f'rewriting contents of folder: {out_path}')
                    warning = False
                if not keep:
                    print(f'{outFile} -> {out_path}{outFile}')
                    os.rename(outFile, out_path + outFile)
            else:
                os.mkdir(out_path)
                if warning:
                    print(f'\ncreating new folder: {out_path}')
                    warning = False
                if not keep:
                    print(f'{outFile} -> {out_path}{outFile}')
                    os.rename(outFile, out_path + outFile)
        print('\n##### CLEANUP DONE #####\n')

    def write_parameters(self):
        keys = self.parameters

        with open('parameters.out', 'w+') as out:
            for key in keys:
                out.write(f'{key}: {self.parameters[key]}\n')

    def run(self):
        from ase.io import read
        from ase.io import write
        from SimLab.simulations import start_sim
        from SimLab.simulations import start_view

        if self.verbosity > 1:
            curr_ase_dftb_command = os.environ["ASE_DFTB_COMMAND"]
            os.environ["ASE_DFTB_COMMAND"] = "dftb+ | tee PREFIX.out"

        molecules = self.fetch_molecule_list()
        molecules.sort()
        for molecule in molecules:
            mol_name = os.path.splitext(os.path.basename(molecule))[0]
            mol = read(molecule)
            out_path = f'{self.parameters["simulation"]}_{self.parameters["method"]}_{mol_name}' + os.sep

            print('Simulation Start')

            if self.parameters["view_only"]:
                start_view(mol, mol_name, out_path, **self.parameters)
            else:
                start_sim(mol, mol_name, out_path, **self.parameters)
                self.write_parameters()
                self.clean_files(out_path)
                start_view(mol, mol_name, out_path, **self.parameters)

        if self.voice:
            os.system('spd-say "It is done."')

        if self.verbosity > 1:
            os.environ["ASE_DFTB_COMMAND"] = curr_ase_dftb_command