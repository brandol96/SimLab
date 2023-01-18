# calculator definitions
def boolean_to_string(parameter):
    # convert logical variable into DFTB+ pattern
    if type(parameter == bool):
        if parameter:
            return 'Yes'
        else:
            return 'No'


def get_grid_origin(mol, n_points):
    A_to_Hr = 0.188972598857892E+01
    positions = mol.get_positions()
    X = []
    Y = []
    Z = []
    for v in positions:
        X.append(v[0])
        Y.append(v[1])
        Z.append(v[2])
    grid_O = [min(X) * A_to_Hr, min(Y) * A_to_Hr, min(Z) * A_to_Hr]
    grid_S = [(max(X) - min(X)) * A_to_Hr / n_points,
              (max(Y) - min(Y)) * A_to_Hr / n_points,
              (max(Z) - min(Z)) * A_to_Hr / n_points]
    return grid_O, grid_S


def fetch_dftb_calc(mol,cluster, **kwargs):
    from ase.calculators.dftb import Dftb
    label = kwargs.get('label')
    kpts = kwargs.get('kpts')
    max_force = kwargs.get('max_force')
    max_driver_steps = kwargs.get('max_driver_steps')
    lattice_opt = boolean_to_string(kwargs.get('lattice_opt'))
    fix_angles = boolean_to_string(kwargs.get('fix_angles'))
    fix_lengths = kwargs.get('fix_lengths').copy()
    n_points = kwargs.get('n_points')
    for i in range(3):
        fix_lengths[i] = boolean_to_string(fix_lengths[i])

    SCC = boolean_to_string(kwargs.get('SCC'))
    max_SCC = kwargs.get('max_SCC')
    max_SCC_steps = kwargs.get('max_SCC_steps')
    fermi_filling = kwargs.get('fermi_filling')
    use_LennardJones = kwargs.get('use_LennardJones')
    grid_O, grid_S = get_grid_origin(mol, n_points)

    eVA_to_HaBohr = 0.01944689673

    if cluster:
        if use_LennardJones:
            calc = Dftb(label=label,
                        Driver_="ConjugateGradient",
                        Driver_MaxForceComponent=max_force * eVA_to_HaBohr,
                        Driver_MaxSteps=max_driver_steps,
                        Driver_MovedAtoms='1:-1',
                        Driver_AppendGeometries='Yes',
                        Hamiltonian_SCC=SCC,
                        Hamiltonian_SCCTolerance=max_SCC,
                        Hamiltonian_MaxSCCIterations=max_SCC_steps,
                        Hamiltonian_Filling=f"Fermi{{Temperature [K] = {fermi_filling} }}",
                        Hamiltonian_Dispersion='LennardJones{Parameters = UFFParameters{}}',
                        Analysis_='',
                        Analysis_ElectrostaticPotential_='',
                        Analysis_ElectrostaticPotential_OutputFile=f'potential_optimize.out',
                        Analysis_ElectrostaticPotential_AppendFile='Yes',
                        Analysis_ElectrostaticPotential_Softening='0.00001',
                        Analysis_ElectrostaticPotential_Grid_='',
                        Analysis_ElectrostaticPotential_Grid_Spacing=f'{grid_S[0]} {grid_S[1]} {grid_S[1]}',
                        Analysis_ElectrostaticPotential_Grid_Origin=f'{grid_O[0]} {grid_O[1]} {grid_O[2]}',
                        Analysis_ElectrostaticPotential_Grid_GridPoints=f'{n_points} {n_points} 1',
                        Analysis_ElectrostaticPotential_Grid_Directions='1 0 0 0 1 0 0 0 1',
                        )
        else:
            calc = Dftb(label=label,
                        Driver_="ConjugateGradient",
                        Driver_MaxForceComponent=max_force * eVA_to_HaBohr,
                        Driver_MaxSteps=max_driver_steps,
                        Driver_MovedAtoms='1:-1',
                        Driver_AppendGeometries='Yes',
                        Hamiltonian_SCC=SCC,
                        Hamiltonian_SCCTolerance=max_SCC,
                        Hamiltonian_MaxSCCIterations=max_SCC_steps,
                        Hamiltonian_Filling=f"Fermi{{Temperature [K] = {fermi_filling} }}",
                        Analysis_='',
                        Analysis_ElectrostaticPotential_='',
                        Analysis_ElectrostaticPotential_OutputFile=f'potential_optimize.out',
                        Analysis_ElectrostaticPotential_AppendFile='Yes',
                        Analysis_ElectrostaticPotential_Softening='0.00001',
                        Analysis_ElectrostaticPotential_Grid_='',
                        Analysis_ElectrostaticPotential_Grid_Spacing=f'{grid_S[0]} {grid_S[1]} {grid_S[1]}',
                        Analysis_ElectrostaticPotential_Grid_Origin=f'{grid_O[0]} {grid_O[1]} {grid_O[2]}',
                        Analysis_ElectrostaticPotential_Grid_GridPoints=f'{n_points} {n_points} 1',
                        Analysis_ElectrostaticPotential_Grid_Directions='1 0 0 0 1 0 0 0 1',
                        )
    else:
        if use_LennardJones:
            calc = Dftb(label=label,
                        kpts=kpts,
                        Driver_="ConjugateGradient",
                        Driver_MaxForceComponent=max_force * eVA_to_HaBohr,
                        Driver_MovedAtoms='1:-1',
                        Driver_LatticeOpt=lattice_opt,
                        Driver_FixAngles=fix_angles,
                        Driver_FixLengths=fix_lengths[0] + fix_lengths[1] + fix_lengths[2],
                        Driver_MaxSteps=max_driver_steps,
                        Driver_AppendGeometries='Yes',
                        Hamiltonian_SCC=SCC,
                        Hamiltonian_SCCTolerance=max_SCC,
                        Hamiltonian_MaxSCCIterations=max_SCC_steps,
                        Hamiltonian_Filling=f"Fermi{{Temperature [K] = {fermi_filling} }}",
                        Hamiltonian_Dispersion='LennardJones{Parameters = UFFParameters{}}',
                        Analysis_='',
                        Analysis_ElectrostaticPotential_='',
                        Analysis_ElectrostaticPotential_OutputFile=f'potential_optimize.out',
                        Analysis_ElectrostaticPotential_AppendFile='Yes',
                        Analysis_ElectrostaticPotential_Softening='0.00001',
                        Analysis_ElectrostaticPotential_Grid_='',
                        Analysis_ElectrostaticPotential_Grid_Spacing=f'{grid_S[0]} {grid_S[1]} {grid_S[1]}',
                        Analysis_ElectrostaticPotential_Grid_Origin=f'{grid_O[0]} {grid_O[1]} {grid_O[2]}',
                        Analysis_ElectrostaticPotential_Grid_GridPoints=f'{n_points} {n_points} 1',
                        Analysis_ElectrostaticPotential_Grid_Directions='1 0 0 0 1 0 0 0 1',
                        )
        else:
            calc = Dftb(label=label,
                        kpts=kpts,
                        Driver_="ConjugateGradient",
                        Driver_MaxForceComponent=max_force * eVA_to_HaBohr,
                        Driver_MovedAtoms='1:-1',
                        Driver_LatticeOpt=lattice_opt,
                        Driver_FixAngles=fix_angles,
                        Driver_FixLengths=self.fixLengths[0] + self.fixLengths[1] + self.fixLengths[2],
                        Driver_MaxSteps=max_driver_steps,
                        Driver_AppendGeometries='Yes',
                        Hamiltonian_SCC=SCC,
                        Hamiltonian_SCCTolerance=max_SCC,
                        Hamiltonian_MaxSCCIterations=max_SCC_steps,
                        Hamiltonian_Filling=f"Fermi{{Temperature [K] = {fermi_filling} }}",
                        Analysis_='',
                        Analysis_ElectrostaticPotential_='',
                        Analysis_ElectrostaticPotential_OutputFile=f'potential_optimize.out',
                        Analysis_ElectrostaticPotential_AppendFile='Yes',
                        Analysis_ElectrostaticPotential_Softening='0.00001',
                        Analysis_ElectrostaticPotential_Grid_='',
                        Analysis_ElectrostaticPotential_Grid_Spacing=f'{grid_S[0]} {grid_S[1]} {grid_S[1]}',
                        Analysis_ElectrostaticPotential_Grid_Origin=f'{grid_O[0]} {grid_O[1]} {grid_O[2]}',
                        Analysis_ElectrostaticPotential_Grid_GridPoints=f'{n_points} {n_points} 1',
                        Analysis_ElectrostaticPotential_Grid_Directions='1 0 0 0 1 0 0 0 1',
                        )
    return calc
