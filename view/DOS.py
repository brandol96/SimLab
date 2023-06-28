import os
import matplotlib.pyplot as plt
from SimLab.utils import read_dos_dftb
from SimLab.utils import read_fermi_levels_dftb

# TODO: make this good
def run_dftb(mol, mol_name, out_path,
             interactive_plot):
    os.system(f'dp_dos {out_path}band.out {out_path}{mol_name}.dos.dat')
    ene, dos = read_dos_dftb(out_path, mol_name)
    homo, lumo, gap, fermi_e = read_fermi_levels_dftb(out_path, mol_name)
    plt.plot(ene, dos)

    if interactive_plot:
        pass
        plt.show()
        #plt.clf()
    else:
        plt.clf()