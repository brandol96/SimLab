import os
import matplotlib.pyplot as plt
from SimLab.utils import read_dos_dftb
from SimLab.utils import read_fermi_levels_dftb


def center_to_fermi(energy, dos, fermi_e, dos_range ,eigen = []):
    # set DOS center to zero
    ene = [Ene - fermi_e for Ene in energy]
    if eigen:
        print('moving eigen')
        eigen = [Lam - fermi_e for Lam in eigen]

    # remove from ene and eigen all values outside dos_range
    aux_dos = []
    aux_ene = []

    for E, D in zip(ene, dos):
        if -dos_range <= E <= dos_range:
            aux_dos.append(D)
            aux_ene.append(E)


    if eigen:
        return aux_ene, aux_dos, eigen
    else:
        return aux_ene, aux_dos


def run_dftb(mol, mol_name, out_path, interactive_plot, dos_range, plot_PDOS = False):
    os.system(f'dp_dos {out_path}band.out {out_path}{mol_name}.dos.dat')

    dos_file = f'{mol_name}.dos.dat'
    ene, dos, eigen, eigen_occs = read_dos_dftb(out_path,dos_file, True)
    homo, lumo, gap, fermi_e = read_fermi_levels_dftb(out_path, mol_name)

    # setup figure
    fig, ax = plt.subplots(figsize=(20, 5))
    fig.suptitle(mol_name.replace("_", " ") + ' DOS', fontsize=20)

    ene, dos, eigen = center_to_fermi(ene, dos, fermi_e, dos_range, eigen)

    # plot stuff
    ax.plot(ene, dos)
    ax.plot([0,0],[min(dos),max(dos)],'--')
    for e in eigen:
        if -dos_range <= e <= dos_range:
            ax.plot([e, e], [0, 0.1 * max(dos)], color='black')

    # text labels with relevant info
    ax.text(min(ene), max(dos), f'Fermi Energy: {fermi_e} eV, gap: {gap} eV', fontsize=12)
    ax.set_ylabel('DOS', fontsize=20)
    ax.set_xlabel('Energy [eV]', fontsize=20)
    ax.set_yticks([])

    shell_dict = {'1':'s',
                  '2':'p',
                  '3':'d',
                  '4':'f',}

    if plot_PDOS:
        out_files = os.listdir(out_path)
        for output in out_files:
            if 'PDOS' in output:
                data = output.split('_')
                data = data[1].split('.')
                atom = data[0]
                shell = data[1]
                orbital = data[2]
                dos_file = f'{out_path}{output}.dos.dat'
                os.system(f'dp_dos {out_path}{output} {dos_file}')
                ene, dos = read_dos_dftb(out_path, dos_file)
                ene, dos = center_to_fermi(ene, dos, fermi_e, dos_range)
                ax.plot(ene, dos, label=f'pdos: {atom} {shell_dict[shell]} {orbital}')

    plt.tight_layout()
    ax.legend()
    fig.savefig(f'{out_path}DFTB_{mol_name}_DOS.png')
    if interactive_plot:
        plt.show()
    else:
        plt.close('all')
