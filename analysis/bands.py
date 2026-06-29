def run(OMP_threads, MPI_cores, verbosity,
        mol, mol_name, sampling):
    import os
    from ase.calculators.dftb import Dftb
    from SimLab.calculator import set_parallelism
    bands = Dftb(atoms=mol,
                 Hamiltonian_KPointsAndWeights=f'KLines {{ \n {sampling} }}',
                 Hamiltonian_SCC='Yes',
                 Hamiltonian_ReadInitialCharges='Yes',
                 Hamiltonian_MaxSCCIterations='100',
                 #Analysis_='',
                 #Analysis_WriteEigenVectors="Yes",
                 #Options_='',
                 #Options_WriteDetailedXml='Yes'
                 )
    # run optimization through DFTB+ implemented routines
    bands = set_parallelism(bands, OMP_threads, MPI_cores, verbosity)
    bands.calculate(mol)

    # dp_band to transform band.out into plottable data
    os.system(f'dp_bands band.out {mol_name}.band')
    os.system(f'dp_dos band.out {mol_name}.dos.dat')
