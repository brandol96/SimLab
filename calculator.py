# calculator definitions
from pydoc import replace


def set_parallelism(calc,OMP_threads,MPI_cores,verbosity):
    import os
    if MPI_cores != 1:
        os.environ["OMP_NUM_THREADS"] = "1"
        os.environ["OPENBLAS_NUM_THREADS"] = "1"
        dftb_mpi_bin = os.environ["DFTB_MPI_BIN"]
        dftb_mpi_lib = os.environ["DFTB_MPI_LIB"]
        inline_env = f"LD_LIBRARY_PATH={dftb_mpi_lib}"
        if verbosity > 2:
            calc.command = f'{inline_env} mpiexec -np {MPI_cores} -x LD_LIBRARY_PATH {dftb_mpi_bin}dftb+ | tee PREFIX.out'
        else:
            calc.command = f'{inline_env} mpiexec -np {MPI_cores} -x LD_LIBRARY_PATH {dftb_mpi_bin}dftb+ > PREFIX.out'
        return calc
    else:
        os.environ["ASE_DFTB_COMMAND"] = f'dftb_omp > PREFIX.out'
        os.environ["OMP_NUM_THREADS"] = str(OMP_threads)
        os.environ["OMP_PROC_BIND"] = 'true'
        os.environ["OMP_PLACES"] = 'cores'
        os.environ["OPENBLAS_NUM_THREADS"] = str(OMP_threads)
        dftb_omp_bin = os.environ["DFTB_OMP_BIN"]
        dftb_omp_lib = os.environ["DFTB_OMP_LIB"]
        inline_env = f"LD_LIBRARY_PATH={dftb_omp_lib}"
        calc.set(Parallel_='',
                 Parallel_UseOmpThreads='Yes')
        print(dftb_omp_bin)
        print(dftb_omp_lib)
        if verbosity > 2:
            calc.command = f'{inline_env} {dftb_omp_bin}dftb+ | tee PREFIX.out'
        else:
            calc.command = f'{inline_env} {dftb_omp_bin}dftb+ > PREFIX.out'

        return calc

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

def fetch_spin_constants(spin_file_path):
    import os
    spin_constants_dict = dict()
    spin_constants_list = []
    dict_key = False
    with open(f'{spin_file_path}{os.sep}spinw.txt','r') as spin_file:
        for line in spin_file:
            line = line.strip()
            if line:
                if ':' in line:
                    if dict_key:
                        spin_constants_dict[dict_key] = spin_constants_list
                        spin_constants_list = []
                    dict_key = line.replace(':','')
                else:
                    spin_constants_list.append(line.split('    '))
    spin_constants_dict[dict_key] = spin_constants_list
    return(spin_constants_dict)


def fetch_dftb_calc(mol, cluster, **kwargs):
    from ase.calculators.dftb import Dftb
    label = kwargs.get('label')
    kpts = kwargs.get('kpts')
    SKFiles = kwargs.get('SKFiles')
    max_force = kwargs.get('max_force')
    max_driver_steps = kwargs.get('max_driver_steps')
    lattice_opt = boolean_to_string(kwargs.get('lattice_opt'))
    fix_angles = boolean_to_string(kwargs.get('fix_angles'))
    spin_polarisation = kwargs.get('spin_polarisation')
    fix_lengths = kwargs.get('fix_lengths').copy()
    n_points = kwargs.get('n_points')
    max_atom_step = kwargs.get('max_atom_step')
    for i in range(3):
        fix_lengths[i] = boolean_to_string(fix_lengths[i])

    SCC = boolean_to_string(kwargs.get('SCC'))
    max_SCC = kwargs.get('max_SCC')
    max_SCC_steps = kwargs.get('max_SCC_steps')
    fermi_filling = kwargs.get('fermi_filling')
    use_LennardJones = kwargs.get('use_LennardJones')
    write_eigens_bin = boolean_to_string(kwargs.get('write_eigens_bin'))
    write_detail_xml = boolean_to_string(kwargs.get('write_detail_xml'))

    #grid_O, grid_S = get_grid_origin(mol, n_points)

    eVA_to_HaBohr = 0.01944689673

    chemSymbs = list(set(mol.get_chemical_symbols()))
    PDOS_string = f'{{\n'
    for chemSymb in chemSymbs:
        PDOS_string += f'Region = {{\n Atoms = {chemSymb}\n OrbitalResolved = Yes\n Label = PDOS_{chemSymb} }} \n'
    PDOS_string += '}'

    # I have to think better about the usage of charges in '.dat' or '.bin'
    # There are cases where each is important

    calc_dict = dict(label=label,
                     Driver_="GeometryOptimisation",
                     Driver_GeometryOptimisation="Lbfgs{ Memory = 20 }",
                     Driver_MaxSteps=max_driver_steps,
                     Driver_Convergence=f"{{ GradMax = {max_force * eVA_to_HaBohr} }}",
                     Driver_MaxAtomStep=max_atom_step,
                     Driver_MovedAtoms='1:-1',
                     Driver_AppendGeometries='Yes',
                     Hamiltonian_SCC=SCC,
                     Hamiltonian_SCCTolerance=max_SCC,
                     Hamiltonian_MaxSCCIterations=max_SCC_steps,
                     Hamiltonian_ReadInitialCharges='No',
                     Hamiltonian_Filling=f"Fermi{{Temperature [K] = {fermi_filling} }}",
                     Analysis_='',
                     Analysis_WriteEigenvectors=write_eigens_bin,
                     Options_WriteChargesAsText='No',
                     Options_WriteDetailedXml=write_detail_xml,
                     ParserOptions_="",
                     ParserOptions_ParserVersion=11,
                     )

    if not cluster:
        calc_dict['kpts'] = kpts

    if use_LennardJones:
        calc_dict['Hamiltonian_Dispersion']='LennardJones{Parameters = UFFParameters{}}'

    if spin_polarisation:
        spin_constants = fetch_spin_constants(SKFiles)
        print('Spin constants enabled!')
        spin_constants_string = '{\n      ShellResolvedSpin = Yes\n'
        for chemSymb in chemSymbs:
            spin_constants_string += f'      {chemSymb} {{'
            for orbital in spin_constants[chemSymb]:
                spin_constants_string += f'{' '.join(orbital)} '
            spin_constants_string += '}\n'
        spin_constants_string += '    }'

        calc_dict['Hamiltonian_SpinConstants'] = spin_constants_string
        calc_dict['Hamiltonian_ShellResolvedSCC'] = 'No'
        calc_dict['Hamiltonian_SpinPolarisation']='Colinear{}'


    calc = Dftb (**calc_dict)

    return calc


def fetch_dftb_calc_deprecated(mol, cluster, **kwargs):
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
                        Hamiltonian_ReadInitialCharges='Yes',
                        Analysis_='',
                        Analysis_WriteEigenvectors='Yes',
                        #Analysis_ElectrostaticPotential_='',
                        #Analysis_ElectrostaticPotential_OutputFile=f'potential_optimize.out',
                        #Analysis_ElectrostaticPotential_AppendFile='Yes',
                        #Analysis_ElectrostaticPotential_Softening='0.01',
                        #Analysis_ElectrostaticPotential_Grid_='',
                        #Analysis_ElectrostaticPotential_Grid_Spacing=f'{grid_S[0]} {grid_S[1]} {grid_S[1]}',
                        #Analysis_ElectrostaticPotential_Grid_Origin=f'{grid_O[0]} {grid_O[1]} {grid_O[2]}',
                        #Analysis_ElectrostaticPotential_Grid_GridPoints=f'{n_points} {n_points} 1',
                        #Analysis_ElectrostaticPotential_Grid_Directions='1 0 0 0 1 0 0 0 1',
                        Options_='',
                        Options_ReadChargesAsText='Yes',
                        Options_WriteDetailedXml='Yes',
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
                        Hamiltonian_ReadInitialCharges='Yes',
                        Analysis_='',
                        Analysis_WriteEigenvectors='Yes',
                        #Analysis_ElectrostaticPotential_='',
                        #Analysis_ElectrostaticPotential_OutputFile=f'potential_optimize.out',
                        #Analysis_ElectrostaticPotential_AppendFile='Yes',
                        #Analysis_ElectrostaticPotential_Softening='0.01',
                        #Analysis_ElectrostaticPotential_Grid_='',
                        #Analysis_ElectrostaticPotential_Grid_Spacing=f'{grid_S[0]} {grid_S[1]} {grid_S[1]}',
                        #Analysis_ElectrostaticPotential_Grid_Origin=f'{grid_O[0]} {grid_O[1]} {grid_O[2]}',
                        #Analysis_ElectrostaticPotential_Grid_GridPoints=f'{n_points} {n_points} 1',
                        #Analysis_ElectrostaticPotential_Grid_Directions='1 0 0 0 1 0 0 0 1',
                        Options_='',
                        Options_ReadChargesAsText='Yes',
                        Options_WriteDetailedXml='Yes',
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
                        Driver_FixLengths=f'{fix_lengths[0]} {fix_lengths[1]} {fix_lengths[2]}',
                        Driver_MaxSteps=max_driver_steps,
                        Driver_AppendGeometries='Yes',
                        Hamiltonian_SCC=SCC,
                        Hamiltonian_SCCTolerance=max_SCC,
                        Hamiltonian_MaxSCCIterations=max_SCC_steps,
                        Hamiltonian_Filling=f"Fermi{{Temperature [K] = {fermi_filling} }}",
                        Hamiltonian_Dispersion='LennardJones{Parameters = UFFParameters{}}',
                        Hamiltonian_ReadInitialCharges='Yes',
                        Analysis_='',
                        Analysis_WriteEigenvectors='Yes',
                        #Analysis_ElectrostaticPotential_='',
                        #Analysis_ElectrostaticPotential_OutputFile=f'potential_optimize.out',
                        #Analysis_ElectrostaticPotential_AppendFile='Yes',
                        #Analysis_ElectrostaticPotential_Softening='0.01',
                        #Analysis_ElectrostaticPotential_Grid_='',
                        #Analysis_ElectrostaticPotential_Grid_Spacing=f'{grid_S[0]} {grid_S[1]} {grid_S[1]}',
                        #Analysis_ElectrostaticPotential_Grid_Origin=f'{grid_O[0]} {grid_O[1]} {grid_O[2]}',
                        #Analysis_ElectrostaticPotential_Grid_GridPoints=f'{n_points} {n_points} 1',
                        #Analysis_ElectrostaticPotential_Grid_Directions='1 0 0 0 1 0 0 0 1',
                        Options_='',
                        Options_ReadChargesAsText='Yes',
                        Options_WriteDetailedXml='Yes',
                        )
        else:
            calc = Dftb(label=label,
                        kpts=kpts,
                        Driver_="ConjugateGradient",
                        Driver_MaxForceComponent=max_force * eVA_to_HaBohr,
                        Driver_MovedAtoms='1:-1',
                        Driver_LatticeOpt=lattice_opt,
                        Driver_FixAngles=fix_angles,
                        Driver_FixLengths=f'{fix_lengths[0]} {fix_lengths[1]} {fix_lengths[2]}',
                        Driver_MaxSteps=max_driver_steps,
                        Driver_AppendGeometries='Yes',
                        Hamiltonian_SCC=SCC,
                        Hamiltonian_SCCTolerance=max_SCC,
                        Hamiltonian_MaxSCCIterations=max_SCC_steps,
                        Hamiltonian_ReadInitialCharges='Yes',
                        Hamiltonian_Filling=f"Fermi{{Temperature [K] = {fermi_filling} }}",
                        Analysis_='',
                        Analysis_WriteEigenvectors='Yes',
                        #Analysis_ElectrostaticPotential_='',
                        #Analysis_ElectrostaticPotential_OutputFile=f'potential_optimize.out',
                        #Analysis_ElectrostaticPotential_AppendFile='Yes',
                        #Analysis_ElectrostaticPotential_Softening='0.01',
                        #Analysis_ElectrostaticPotential_Grid_='',
                        #Analysis_ElectrostaticPotential_Grid_Spacing=f'{grid_S[0]} {grid_S[1]} {grid_S[1]}',
                        #Analysis_ElectrostaticPotential_Grid_Origin=f'{grid_O[0]} {grid_O[1]} {grid_O[2]}',
                        #Analysis_ElectrostaticPotential_Grid_GridPoints=f'{n_points} {n_points} 1',
                        #Analysis_ElectrostaticPotential_Grid_Directions='1 0 0 0 1 0 0 0 1',
                        Options_='',
                        Options_ReadChargesAsText='Yes',
                        Options_WriteDetailedXml='Yes',
                        )
    return calc
