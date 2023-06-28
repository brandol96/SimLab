import numpy as np
import matplotlib.pyplot as plt
from SimLab.utils import path_dftb
from SimLab.utils import read_dos_dftb
from SimLab.utils import read_fermi_levels_dftb
from SimLab.utils import output_fermi_levels_dftb

def run_dftb(mol, mol_name, out_path, zoom, path, dK,
             interactive_plot):
    band_data = np.genfromtxt(f'{out_path}{mol_name}.band_tot.dat')
    # band_data contains nXY plot of the band structure, it is a matrix of the form
    # [[kpt1 E_band1 E_band2 ... E_bandM]
    #  [kpt2 E_band1 E_band2 ... E_bandM]
    #  ...
    #  [kptN E_band1 E_band2 ... E_bandM]]

    ene, dos = read_dos_dftb(out_path, mol_name)
    # read_dos_dftb return two lists with the aligned plot data
    # ene = [Ene0, Ene1, Ene2, ..., Enen]
    # dos = [Dos0, Dos1m Dos2, ..., Dosn]

    homo, lumo, gap, fermi_e = read_fermi_levels_dftb(out_path, mol_name)
    # homo -> info about found homo [KPT,BAND,EV]
    # lumo -> info about found lumo [KPT,BAND,EV]
    # gap  -> gap value calculated for current file
    # fermi_e -> fermi energy calculated for current file

    output_fermi_levels_dftb(mol_name, homo, lumo, fermi_e, gap)
    # output fermi level without replacing old results

    # center plot into fermi energy
    for i in range(len(band_data)):
        band_data[i][1:] = band_data[i][1:] - fermi_e
    for i in range(len(dos)):
        ene[i] = ene[i] - fermi_e
#    fermi_e = 0

    # setup figure
    fig = plt.figure(1, figsize=(8, 10))  # start a figure
    fig.suptitle(mol_name.replace("_", " ") + ' band structure')

    # bands axes
    ax = fig.add_axes([.12, .07, .67, .85])  # axes [left, bottom, width, height]
    ax.set_xticks([])
    ax.set_ylabel('$E - E_f$ (eV)')
    if zoom > 0:
        ax.set_ylim([-zoom, +zoom])

    # dos axes
    dosax = fig.add_axes([.8, .07, .17, .85])  # axes [left, bottom, width, height]
    dosax.fill_between(dos, ene)
    dosax.set_yticks([])
    dosax.set_xticks([])
    dosax.set_xlabel("DOS")
    if zoom > 0:
        dosax.set_ylim([-zoom, +zoom])

    # plot bands and dos
    for i in range(band_data.shape[1] - 1):
        ax.plot(band_data[:, 0], band_data[:, i + 1])
    dosax.plot(dos, ene, color='black')

    # plot vertical lines to indicate path taken

    # find Y-points to plot de KLines of the greater and lower energy values on Y-axis
    # to plot the lines in visual range
    ymin = 0.0
    ymax = 0.0
    for i in range(len(band_data[0, 1:])):
        if min(band_data[:, i + 1]) < ymin:
            ymin = min(band_data[:, i + 1])
        if max(band_data[:, i + 1]) > ymax:
            ymax = max(band_data[:, i + 1])

    dftb_path = path_dftb(path, dK, mol, 0, False)
    data = dftb_path.split('\n')

    # the strategy here:
    # user provided 'path' is aligned to the generated 'dftb_path'
    # I'll just exploit this to plot the lines the first letter in 'path'
    # will allways be the first point in the 'dftb_path' entry and so on...

    path_len = len(data)-1 # my code generates and extra '\n' at the end
    position = 0
    for i in range(path_len):
        point = path[i]
        position += int(data[i].split()[0])
        ax.plot([position, position], [ymin, ymax], color='black')
        if zoom <= 0:
            ax.text(position, ymin, point)
        else:
            ax.text(position, -zoom, point)

    fig.savefig(f'{out_path}DFTB_{mol_name}_Bands.png')
    if interactive_plot:
        plt.show()
    else:
        plt.clf()
