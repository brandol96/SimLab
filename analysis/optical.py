def get_grid_origin(mol,n_points):
    A_to_Hr = 0.188972598857892E+01
    positions = mol.get_positions()
    X = []
    Y = []
    Z = []
    for v in positions:
        X.append(v[0])
        Y.append(v[1])
        Z.append(v[2])
    grid_O = [min(X)*A_to_Hr, min(Y)*A_to_Hr, min(Z)*A_to_Hr]
    grid_S = [(max(X)-min(X))*A_to_Hr/n_points,
              (max(Y)-min(Y))*A_to_Hr/n_points,
              (max(Z)-min(Z))*A_to_Hr/n_points]
    return grid_O, grid_S


def run_kick(mol, max_SCC, max_SCC_steps, fermi_filling,
             total_time, time_step, field_strength, n_points, direction):
    from ase.calculators.dftb import Dftb

    # writing DFTB electron dynamics manually
    totalSteps = int(total_time / time_step)

    electron_dynamics = '{'
    electron_dynamics += f'\nSteps = {totalSteps}'
    electron_dynamics += f'\nTimeStep [fs] = {time_step}'
    electron_dynamics += f'\nPerturbation = Kick {{ \nPolarizationDirection = {direction}\n }}'
    electron_dynamics += f'\nFieldStrength [V/A] = {field_strength}'
    electron_dynamics += f'\nWriteEnergyAndCharges = Yes}}'

    grid_O, grid_S = get_grid_origin(mol, n_points)
    optical = Dftb(atoms=mol,
                   label=f'optical_kick_{direction}',
                   Hamiltonian_SCC='Yes',
                   Hamiltonian_SCCTolerance=max_SCC,
                   Hamiltonian_ReadInitialCharges='Yes',
                   Hamiltonian_MaxSCCIterations=max_SCC_steps,
                   Hamiltonian_Filling=f"Fermi{{Temperature [K] = {fermi_filling} }}",
                   ElectronDynamics=electron_dynamics,
                   Analysis_='',
                   Analysis_ElectrostaticPotential_='',
                   Analysis_ElectrostaticPotential_OutputFile=f'potential_kick_{direction}.out',
                   Analysis_ElectrostaticPotential_AppendFile='Yes',
                   Analysis_ElectrostaticPotential_Softening='0.1',
                   Analysis_ElectrostaticPotential_Grid_='',
                   Analysis_ElectrostaticPotential_Grid_Spacing=f'{grid_S[0]} {grid_S[1]} {grid_S[1]}',
                   Analysis_ElectrostaticPotential_Grid_Origin=f'{grid_O[0]} {grid_O[1]} {grid_O[2]}',
                   Analysis_ElectrostaticPotential_Grid_GridPoints=f'{n_points} {n_points} 1',
                   Analysis_ElectrostaticPotential_Grid_Directions='1 0 0 0 1 0 0 0 1',
                   Options_='',
                   Options_WriteChargesAsText='Yes')

    # run calculation through DFTB+ implemented routines
    optical.calculate(mol)


def run_laser(mol, max_SCC, max_SCC_steps, fermi_filling,
              total_time, time_step, laser_energy, field_strength, n_points, direction):
    from ase.calculators.dftb import Dftb

    # writing DFTB electron dynamics manually
    totalSteps = int(total_time / time_step)
    laser_pol = f'{direction[0]} {direction[1]} {direction[2]}'

    electron_dynamics = '{'
    electron_dynamics += f'\nSteps = {totalSteps}'
    electron_dynamics += f'\nTimeStep [fs] = {time_step}'
    electron_dynamics += f'\nPerturbation = Laser {{ \nPolarizationDirection = {laser_pol}'
    electron_dynamics += f'\nLaserEnergy [eV] ={laser_energy} }}'
    electron_dynamics += f'\nFieldStrength [V/A] = {field_strength}'
    electron_dynamics += f'\nWriteFrequency = 20'
    # electron_dynamics += f'\nIonDynamics = Yes'
    # electron_dynamics += f'\nInitialTemperature [K] = 0.0'
    # electron_dynamics += f'\nPopulations = Yes'
    electron_dynamics += f'\nWriteEnergyAndCharges = Yes}}'

    grid_O, grid_S = get_grid_origin(mol,n_points)
    optical = Dftb(atoms=mol,
                   label=f'optical_laser_{direction[0]}{direction[1]}{direction[2]}',
                   Hamiltonian_SCC='Yes',
                   Hamiltonian_SCCTolerance=max_SCC,
                   Hamiltonian_ReadInitialCharges='Yes',
                   Hamiltonian_MaxSCCIterations=max_SCC_steps,
                   Hamiltonian_Filling=f"Fermi{{Temperature [K] = {fermi_filling} }}",
                   ElectronDynamics=electron_dynamics,
                   Analysis_='',
                   Analysis_ElectrostaticPotential_='',
                   Analysis_ElectrostaticPotential_OutputFile=f'potential_laser_{direction[0]}{direction[1]}{direction[2]}.out',
                   Analysis_ElectrostaticPotential_AppendFile='Yes',
                   Analysis_ElectrostaticPotential_Softening='0.00001',
                   Analysis_ElectrostaticPotential_Grid_='',
                   Analysis_ElectrostaticPotential_Grid_Spacing=f'{grid_S[0]} {grid_S[1]} {grid_S[1]}',
                   Analysis_ElectrostaticPotential_Grid_Origin=f'{grid_O[0]} {grid_O[1]} {grid_O[2]}',
                   Analysis_ElectrostaticPotential_Grid_GridPoints=f'{n_points} {n_points} 1',
                   Analysis_ElectrostaticPotential_Grid_Directions='1 0 0 0 1 0 0 0 1',
                   Options_='',
                   Options_WriteChargesAsText='Yes')

    # run calculation through DFTB+ implemented routines
    optical.calculate(mol)
