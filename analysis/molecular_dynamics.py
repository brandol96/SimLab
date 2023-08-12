def highest_mode():
    search_flag = False
    with open('modes.out') as modes_file:
        for line in modes_file:
            if 'cm-1' in line:
                search_flag = True
            if search_flag:
                data = line.split()
                if len(data)==2:
                    data = line.split()
                    vib_mode = data[1]
        return vib_mode


def run(mol, mol_name, kpts, thermostat, temp_profile, time_step,
        SCC, max_SCC, max_SCC_steps, fermi_filling):
    from ase.calculators.dftb import Dftb

    # if this condition if True, then the variable TempProfile has None
    if not temp_profile:
        print('please provide an adequate profile:\n[["keyword1","steps1","temp1"],'
              '\n["keyword2","steps2","temp2"],'
              '\n["keyword3","steps3","temp3"]]')

    # write temperatureProfile
    TempProfileTxt = 'TemperatureProfile {\n'
    for i in range(len(temp_profile)):
        temp_change = temp_profile[i][0]
        temp_step = int(temp_profile[i][1])
        temp_value = float(temp_profile[i][2])*0.316681534524639E-05
        TempProfileTxt += f'{temp_change}    {temp_step}     {temp_value}\n'
    TempProfileTxt += '\n}\n'

    pbc = mol.get_pbc()
    if True in pbc:
        print('Some direction has pbc, be sure kpts provided are correct!')
        md = Dftb(atoms=mol,
                  kpts=kpts)
    else:
        md = Dftb(atoms=mol)

    # setup DFTB+ Calculator
    if thermostat == 'NVE':
        md.set(Driver_="VelocityVerlet",
               Driver_OutputPrefix='NVE',
               Driver_MovedAtoms='1:-1',
               Driver_MDRestartFrequency='10',
               Driver_timestep=time_step*0.413413733365614E+02,
               Driver_Thermostat_='Berendsen',
               Driver_Thermostat_AdaptFillingTemp='Yes',
               Driver_Thermostat_CouplingStrength='1.0',
               Driver_Thermostat_Temperature=TempProfileTxt,
               Hamiltonian_='DFTB',
               Hamiltonian_SCC='No',
               Hamiltonian_Dispersion='LennardJones{Parameters = UFFParameters{}}',
               )
    if thermostat == 'NVT':
        # TODO: modes calculation can be place independently from MD into simulations.py
        # BTW EVERYTHING BEYOND THIS POINT IS UGLY! PLEASE LOOK AWAY!
        modes = Dftb(atoms=mol,
                     label='hessian_run',
                     Driver_="SecondDerivatives",
                     Driver_Delta='1E-4',
                     Hamiltonian_SCC='Yes',
                     Hamiltonian_SCCTolerance=max_SCC,
                     Hamiltonian_MaxSCCIterations=max_SCC_steps,
                     Hamiltonian_Filling=f"Fermi{{Temperature [K] = {fermi_filling} }}",
                     Hamiltonian_Dispersion='LennardJones{Parameters = UFFParameters{}}',
                     )
        mol.set_calculator(modes)
        modes.calculate(mol)
        # TODO: find a good solution to the problem bellow
        # after the hessian is calculated, another part of DFTB+ is called: modes_in.hsd
        # as far as I know, this is not implemented into ASE's DFTB+ calculator because
        # they implemented their own routines. A quick and terribly dirty solution:
        # input files for Modes' code is fairly simpel and should remain unchanged from
        # problem to problem without great complications, I'll just past it inside SimLab
        # and copy each time I want to run this guy
        import shutil
        import os
        import numpy as np
        print('USING MODES CODE! REMEMBER TO CHANGE SKFILES PATH IN THE UGLY FOLDER!!!!')
        print('or do something about it ffs...')
        source = '/home/rbrandolt/python-libs/SimLab/ugly_solutions/modes_in.hsd'
        destination = os.getcwd()
        shutil.copy(source, destination)
        os.system('modes > modes.out')
        vib_mode_freq = highest_mode()
        print(f'found frequency: {vib_mode_freq} cm-1')

        # bellow is a solution to inserting the desired unit [cm-1]
        # in the declaration of coupling strength
        # please fiz this someday into something at least passable
        TempProfileTxt+=f'CouplingStrength [cm^-1] = {vib_mode_freq}\n'

        md.set(Driver_="VelocityVerlet",
               label='NVT_run.out',
               Driver_OutputPrefix='NVT',
               Driver_MovedAtoms='1:-1',
               Driver_MDRestartFrequency='10',
               Driver_timestep=time_step*0.413413733365614E+02,
               Driver_Thermostat_='NoseHoover',
               Driver_Thermostat_AdaptFillingTemp='Yes',
               Driver_Thermostat_Temperature=TempProfileTxt,
               Hamiltonian_='DFTB',
               Hamiltonian_SCC='No',
               Hamiltonian_Dispersion='LennardJones{Parameters = UFFParameters{}}',
               )

    mol.set_calculator(md)
    md.calculate(mol)
