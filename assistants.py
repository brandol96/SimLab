# This contains the basic class of Cody, used to control all of the performed simulations
# every folder file juggling should be places here to be done and written once.
import os
import shutil

# noinspection PyUnresolvedReferences
class Cody:
    # This method sets required initial variables
    # I'll put the simulation parameters here
    def __init__(self, **kwargs):
        # internal setup parameters
        self.valid_simulations = ['optimize','molecular_dynamics','bands','optical','effMass','orbitals']
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
            dos_range=kwargs.get('dos_range', 100),
            verbosity=kwargs.get('verbosity', 0),

            # DFTB stuff
            max_driver_steps=kwargs.get('max_driver_steps', 10000),
            SCC=kwargs.get('SCC', True),
            max_SCC=kwargs.get('max_SCC', 1E-2),
            max_SCC_steps=kwargs.get('max_SCC_steps', 1000),
            fermi_filling=kwargs.get('fermi_filling', 0.0),
            use_LennardJones=kwargs.get('use_LennardJones', False),
            SKFiles=kwargs.get('SKFiles', 'Please supply a choice!'),
            target_orbirals=kwargs.get('target_orbirals', 0),
            WP_grid=kwargs.get('WP_grid', 10),
            WP_Box_View=kwargs.get('WP_Box_View', [1,1,1]),

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
    def fetch_folders_list(self):
        input_list = os.listdir()
        folder_list = []

        for inputFile in input_list:
            if os.path.isdir(inputFile):
                folderParts = inputFile.split('_')
                if folderParts[0] in self.valid_simulations:
                    folder_list.append(inputFile)
        folder_list.sort()
        return input_list, folder_list

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
                   'effMass.out' == outFile or \
                   'effMass.csv' == outFile or \
                   '.vmd' in outFile

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

    def compile_results(self):
        from ase.io import read, write
        from SimLab.utils import read_fermi_levels_dftb
        from SimLab.utils import read_detailed_dftb
        from SimLab.utils import silly_method_to_get_plane_distance
        import csv


        # check for results folder
        all_folders, sim_folders = self.fetch_folders_list(self)
        compiledResults_folder = 'compiled_results'
        if compiledResults_folder in all_folders:
            print(f'found "{compiledResults_folder}" folder')
        else:
            print(f'creating "{compiledResults_folder}" folder')
            os.mkdir(compiledResults_folder)
            os.mkdir(compiledResults_folder + os.sep + 'band')
            os.mkdir(compiledResults_folder + os.sep + 'geom')
            os.mkdir(compiledResults_folder + os.sep + 'traj')

        # startup csv file
        csv_filePath = compiledResults_folder + os.sep + 'compiled_results.csv'
        csv_file = open(csv_filePath, 'w+')

        fieldnames = ['Molecule',
                      'Cell Vector Lenght 1 (Angs)',
                      'Cell Vector Lenght 2 (Angs)',
                      'Cell Vector Lenght 3 (Angs)',
                      'Cell Vector Angle 1 (deg)',
                      'Cell Vector Angle 2 (deg)',
                      'Cell Vector Angle 3 (deg)',
                      'Key Distance 1 (Angs)',
                      'Key Distance 2 (Angs)',
                      'Key Distance 3 (Angs)',
                      'Key Distance 4 (Angs)',
                      'Energy (eV)',
                      'Fermi Level (eV)',
                      'Energy Gap (eV)',
                      'HOMO Band Index',
                      'LUMO Band Index',
                      'Electron Effective Mass',
                      'Hole Effective Mass']
        csv_dict = {'Molecule':'name',
                    'Cell Vector Lenght 1 (Angs)':0.0,
                    'Cell Vector Lenght 2 (Angs)':0.0,
                    'Cell Vector Lenght 3 (Angs)':0.0,
                    'Cell Vector Angle 1 (deg)':0.0,
                    'Cell Vector Angle 2 (deg)':0.0,
                    'Cell Vector Angle 3 (deg)':0.0,
                    'Key Distance 1 (Angs)':0.0,
                    'Key Distance 2 (Angs)':0.0,
                    'Key Distance 3 (Angs)':0.0,
                    'Key Distance 4 (Angs)':0.0,
                    'Energy (eV)':0.0,
                    'Fermi Level (eV)':0.0,
                    'Energy Gap (eV)':0.0,
                    'HOMO Band Index':0.0,
                    'LUMO Band Index':0.0,
                    'Electron Effective Mass':0.0,
                    'Hole Effective Mass':0.0}
        csv_line = []

        csv_writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        csv_writer.writeheader()

        # go trough each valid simulation type and extract the results
        # represented in image files
        for folder in sim_folders:
            if folder.split('_')[0] in self.valid_simulations:
                resultsFiles = os.listdir(folder)
                for fileName in resultsFiles:
                    # very basic test to get .png images
                    if '.png' in fileName:
                        if 'Bands' in fileName or 'EffMass' in fileName:
                            filePath = folder+os.sep+fileName
                            copyPath = (compiledResults_folder + os.sep +
                                        'bands' + os.sep + fileName)
                            shutil.copyfile(filePath, copyPath)
                        if 'end' in fileName:
                            filePath = folder+os.sep+fileName
                            copyPath = (compiledResults_folder + os.sep +
                                        'geom' + os.sep + fileName)
                            shutil.copyfile(filePath, copyPath)
                    if '.traj' in fileName:
                        # make a copy of the dict to avoid updating everything
                        # due to python's variable linking
                        csv_dict_temp = csv_dict.copy()

                        # get molecule name and clean it up a bit
                        molName = os.path.splitext(os.path.basename(fileName))[0]
                        molName = molName.replace('_end','')
                        csv_dict_temp['Molecule'] = molName

                        # read molecule and get its cell
                        end_mol = read(folder+os.sep+fileName)
                        end_cell = end_mol.get_cell_lengths_and_angles()
                        write(f'{mol_name}.traj', end_mol) # write trajectory before expanding
                        end_mol *= (2,2,1) # use 4 unit cells to ensure there are at least 4 atoms
                        csv_dict_temp['Cell Vector Lenght 1 (Angs)'] = round(end_cell[0],3)
                        csv_dict_temp['Cell Vector Lenght 2 (Angs)'] = round(end_cell[1],3)
                        csv_dict_temp['Cell Vector Lenght 3 (Angs)'] = round(end_cell[2],3)
                        csv_dict_temp['Cell Vector Angle 1 (deg)'] = round(end_cell[3],2)
                        csv_dict_temp['Cell Vector Angle 2 (deg)'] = round(end_cell[4],2)
                        csv_dict_temp['Cell Vector Angle 3 (deg)'] = round(end_cell[5],2)

                        # TODO: find out gow to make this an input instead of having to edit the source code
                        # TODO: also, bond lenght being dependent on atom index is very sensible to
                        # TODO geometry construction, it should be more generic so I dot keep this mind
                        csv_dict_temp['Key Distance 1 (Angs)'] = round(end_mol.get_distance(0, 1),3)
                        csv_dict_temp['Key Distance 2 (Angs)'] = round(end_mol.get_distance(0, 2),3)
                        csv_dict_temp['Key Distance 3 (Angs)'] = round(end_mol.get_distance(0, 3),3)
                        print(folder)
                        csv_dict_temp['Key Distance 4 (Angs)'] = silly_method_to_get_plane_distance(end_mol)

                        # get eletronic data
                        bandData_path = folder + os.sep
                        Total_Energy = read_detailed_dftb(folder+os.sep+'detailed.out')
                        homo, lumo, gap, fermi_e = read_fermi_levels_dftb(bandData_path, molName, self.verbosity)
                        csv_dict_temp['Energy (eV)'] = Total_Energy
                        csv_dict_temp['Fermi Level (eV)'] = fermi_e
                        csv_dict_temp['Energy Gap (eV)'] = gap
                        csv_dict_temp['HOMO Band Index'] = homo[1]
                        csv_dict_temp['LUMO Band Index'] = lumo[1]

                        # read effective mass from csv output
                        effMassData_path = '_'.join(['bands'] + folder.split('_')[1:])
                        try:
                            with open(effMassData_path+os.sep+'effMass.csv') as effMass_File:
                                lines = effMass_File.readlines()
                                if lines:  # Check if the file is not empty
                                    last_line = lines[-1].strip()  # .strip() removes trailing newline characters
                                    last_line = last_line.split(',')
                                    csv_dict_temp['Hole Effective Mass'] = last_line[4]
                                    csv_dict_temp['Electron Effective Mass'] = last_line[5]
                                else:
                                    print("File is empty.")
                        except FileNotFoundError:
                            print('Effective mass .csv data not found in Bands folder')



                        csv_line.append(csv_dict_temp)
        csv_writer.writerows(csv_line)
        csv_file.close()


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

        if self.verbosity > 2:
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
                print(os.environ["ASE_DFTB_COMMAND"])
                start_sim(mol, mol_name, out_path, **self.parameters)
                self.write_parameters()
                if self.parameters['simulation'] != 'effMass':
                    self.clean_files(out_path)
                start_view(mol, mol_name, out_path, **self.parameters)

        if self.voice:
            os.system('spd-say "It is done."')

        if self.verbosity > 2:
            os.environ["ASE_DFTB_COMMAND"] = curr_ase_dftb_command
