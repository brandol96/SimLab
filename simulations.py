# here are the set of stuff that a simulation needs done before starting
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
                dftb = fetch_dftb_calc(cluster=False, **kwargs)
            else:
                print('No direction has pbc, DFTB will NOT perform lattice optimization')
                dftb = fetch_dftb_calc(cluster=True, **kwargs)
            # run optimization through DFTB+ implemented routines
            mol.set_calculator(dftb)
            dftb.calculate(mol)

            mol = read(f'geo_end.gen')
            mol.center()
            write(f'{method}_{mol_name}_end.traj', mol)

        else:
            print(f'{simulation} is not possible with {method}')
