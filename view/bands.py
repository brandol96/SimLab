import numpy as np
import matplotlib.pyplot as plt
from SimLab.utils import path_dftb
from SimLab.utils import read_dos_dftb
from SimLab.utils import read_fermi_levels_dftb
from SimLab.utils import output_fermi_levels_dftb


def run_dftb(mol, mol_name, out_path, zoom, path, dK,
             interactive_plot, verbosity):
    fontsize = 30

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

    homo, lumo, gap, fermi_e = read_fermi_levels_dftb(out_path, mol_name, verbosity)
    # homo -> info about found homo [KPT,BAND,EV]
    # lumo -> info about found lumo [KPT,BAND,EV]
    # gap  -> gap value calculated for current file
    # fermi_e -> fermi energy calculated for current file

    # output fermi level without replacing old results
    output_fermi_levels_dftb(mol_name, homo, lumo, fermi_e, gap)

    if lumo[1] == homo[1]:
        print(f'Band {homo[1]} is half-filled! I\'ll use band {homo[1]+1} as LUMO instead!')
        # to select specific band we use the array slicing: band_data[ : , Band_idx ]
        lumo_e = min(band_data[:,lumo[1]+1])
        KPTS, BANDS = np.where(band_data == lumo_e) # find lowest energy KPT
        if verbosity >=3: # print same erngy KPTS if the user desires
            print(f'[KPTS,BANDS] idx found with the same energy: {lumo_e}')
            for kpt, band in zip(KPTS,BANDS):
                print(f'[{kpt},{band}]')

        # take averages to ensure no shenanigans!
        kpt = int(np.average(KPTS))
        band = np.average(BANDS) # OBSERVE: IF THIS RETURNS A NON INTEGER THERE'S SOMETHING AMISS!

        # now update everything
        lumo = [kpt,band,lumo_e]
        fermi_e = round((lumo[2] + homo[2]) / 2, 6)
        gap = round(lumo[2] - homo[2], 6)

        if verbosity >= 2:
            print(f'\n{"molName":<15} {"homo[kpt, Band, eV]":>20} {"lumo[kpt, Band, eV]":>20} {"gap":<3} {"fermi_e":<6}')
            homo_string = f'[{homo[0]},{homo[1]},{homo[2]}]'
            lumo_string = f'[{lumo[0]},{lumo[1]},{lumo[2]}]'
            print(f'{mol_name+"*":<15} {homo_string:>20} {lumo_string:>20} {gap:<3} {fermi_e:<3}')

        # output new results
        output_fermi_levels_dftb(mol_name+'*', homo, lumo, fermi_e, gap)

    # center plot into fermi energy
    for i in range(len(band_data)):
        band_data[i][1:] = band_data[i][1:] - fermi_e
    for i in range(len(dos)):
        ene[i] = ene[i] - fermi_e
    center     = 0.0

    # setup figure
    fig = plt.figure(1, figsize=(10, 12.5))  # start a figure
    #fig.suptitle(mol_name.replace("_", " ") + ' band structure')

    # bands axes
    ax = fig.add_axes([.12, .07, .67, .85])  # axes [left, bottom, width, height]
    ax.set_xticks([])
    ax.set_ylabel(f'$E - E_f$ (eV) |  gap: {gap} eV',fontsize=fontsize)
    if zoom > 0:
        ax.set_ylim([-zoom, +zoom])

    # dos axes
    dosax = fig.add_axes([.8, .07, .17, .85])  # axes [left, bottom, width, height]
    dosax.fill_between(dos, ene, color = 'black')
    dosax.set_yticks([])
    dosax.set_xticks([])
    dosax.set_xlabel("DOS",fontsize=fontsize)
    if zoom > 0:
        dosax.set_ylim([-zoom, +zoom])

    # plot bands and dos
    for i in range(band_data.shape[1] - 1):
        ax.plot(band_data[:, 0], band_data[:, i + 1],color='#272222')
    dosax.plot(dos, ene, color='#272222')

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

    path_len = len(data) - 1  # my code generates and extra '\n' at the end
    position = 0
    for i in range(path_len):
        point = path[i]
        position += int(data[i].split()[0])
        ax.plot([position, position], [ymin, ymax], color='#272222')
        if zoom <= 0:
            ax.text(position, ymin, point, fontsize=fontsize, horizontalalignment='center')
        else:
            ax.text(position, -zoom -0.6, point,fontsize=fontsize, horizontalalignment='center')
    ax.plot([min(band_data[:,0]),max(band_data[:,0])],[center,center], '--', color='#272222')

    # this here improvises a title
    title = mol_name
    ax.text(min(band_data[:, 0]), zoom + 0.15, title, fontsize=fontsize)
    ax.tick_params(axis='both', which='major', labelsize=fontsize)

    fig.savefig(f'{out_path}DFTB_{mol_name}_Bands.png')
    if interactive_plot:
        plt.show()
    else:
        plt.clf()
