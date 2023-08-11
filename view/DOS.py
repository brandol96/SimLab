import os
import matplotlib.pyplot as plt
from SimLab.utils import read_dos_dftb
from SimLab.utils import read_fermi_levels_dftb


def run_dftb(mol, mol_name, out_path, interactive_plot, dos_range):
    os.system(f'dp_dos {out_path}band.out {out_path}{mol_name}.dos.dat')
    ene, dos, eigen = read_dos_dftb(out_path, mol_name, True)
    homo, lumo, gap, fermi_e = read_fermi_levels_dftb(out_path, mol_name)

    # setup figure
    fig, ax = plt.subplots(figsize=(20, 5))
    fig.suptitle(mol_name.replace("_", " ") + ' DOS', fontsize=20)

    # set DOS center to zero
    ene = [Ene - fermi_e for Ene in ene]
    eigen = [Lam - fermi_e for Lam in eigen]

    # remove from ene and eigen all values outside dos_range
    dos_range = [-dos_range, dos_range]
    aux_dos = []
    aux_ene = []

    for E, D in zip(ene, dos):
        if dos_range[0] <= E <= dos_range[1]:
            aux_dos.append(D)
            aux_ene.append(E)
    ene = aux_ene.copy()
    dos = aux_dos.copy()
    del aux_dos, aux_ene

    # plot stuff
    ax.plot(ene, dos)
    for e in eigen:
        if dos_range[0] <= e <= dos_range[1]:
            ax.plot([e, e], [0, 0.1 * max(dos)], color='black')

    # text labels with relevant info
    ax.text(min(ene), max(dos), f'Fermi Energy: {fermi_e} eV, gap: {gap} eV', fontsize=12)
    ax.set_ylabel('DOS', fontsize=20)
    ax.set_xlabel('Energy [eV]', fontsize=20)

    plt.tight_layout()
    fig.savefig(f'{out_path}DFTB_{mol_name}_DOS.png')
    if interactive_plot:
        pass
        plt.show()
        # plt.clf()
    else:
        plt.clf()
