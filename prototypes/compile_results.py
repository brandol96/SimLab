def run():
    # I'll just place here the example code I've made to search bondlenghts based on ASE
    # some oter day I might go back to work this into Cody

    from SimLab.prototypes import DFT_gridGraph
    from ase.geometry.analysis import Analysis
    from ase.visualize import view
    from ase.io import read
    import matplotlib.pyplot as plt
    import numpy as np

    verbose = False
    show_graph = True
    show_mol = True
    step = 0
    X = []
    x = 0
    time = []
    fmax = []
    epot = []
    BFGS_out_fileName = 'hBN_O_Nitrogen_5.100_5.100_30.0_BFGS04'
    with open(f'{BFGS_out_fileName}.log', 'r') as logFile:
        for line in logFile:
            if 'BFGS' in line:
                line = line.split(']')[1]
                char_temp = ''
                for c in line:
                    if c != ' ':
                        char_temp += c
                    if c == ' ':
                        if char_temp:
                            if step == 0:
                                if verbose: print('step start: ' + char_temp)
                            if step == 1:
                                if verbose: print('step:  ' + str(x))
                                X.append(x)
                                x += 1

                                if verbose: print('Energy: ' + char_temp)
                                epot.append(float(char_temp))
                            del char_temp
                            char_temp = ''
                            step += 1
                step = 0
                if verbose: print('fmax: ' + char_temp)
                fmax.append(float(char_temp))

    min_epot_idx = epot.index(min(epot))
    min_fmax_idx = fmax.index(min(fmax))
    print(BFGS_out_fileName)
    print('minimal energy step: ', min_epot_idx, epot[min_epot_idx], fmax[min_epot_idx])
    print('minimal force step: ', min_fmax_idx, epot[min_fmax_idx], fmax[min_fmax_idx])

    # print lattice and bonds
    mol = read(f'{BFGS_out_fileName}.traj', min_epot_idx)
    ana = Analysis(mol)

    print(mol.cell.cellpar())
    symbols = mol.symbols
    for sym in symbols:
        if sym != 'B' and sym != 'N':
            Ads = sym
            break
    CCBonds = ana.get_bonds('B', 'N')
    CCBondValues = ana.get_values(CCBonds)
    # print(CCBondValues)

    # get adsorbate distance from sheet and carbon puck distance
    Z_list = []
    search_idxs = np.concatenate((mol.symbols.indices()['B'], mol.symbols.indices()['N']))
    search_idxs = np.sort(search_idxs)
    for BN_idx in search_idxs:
        position = mol[BN_idx].position
        Z_list.append(position[2])

    z_avg = np.average(Z_list)
    z_std = np.std(Z_list)
    print(f'Start Average Z = {z_avg:.3f} ± {z_std:.3f}')
    pop_idx_list = []
    for position in Z_list:
        if verbose: print(Z_list.index(position), position)
        if position > z_avg + z_std / 2 or position < z_avg - z_std / 2:
            idx = Z_list.index(position)
            pop_idx_list.append(idx)

    for index in sorted(pop_idx_list, reverse=True):
        print(f'Atom {index} is pucked out by {abs(Z_list[index] - z_avg):.3f}')
        del Z_list[index]

    z_avg = np.average(Z_list)
    z_std = np.std(Z_list)
    print(f'New Average Z for sheet: {z_avg:.3f} ± {z_std:.3f}')

    # get adsorbate information based on hetereatom presence
    if 'Ads' in globals():
        BAdsBonds = ana.get_bonds('B', Ads)
        if all(sublist for sublist in BAdsBonds):
            BAdsBondValues = ana.get_values(BAdsBonds)
            print(f'B-Ads bond lenght: {np.average(BAdsBondValues):.3f} Angs')

        NAdsBonds = ana.get_bonds('N', Ads)
        if all(sublist for sublist in NAdsBonds):
            NAdsBondValues = ana.get_values(NAdsBonds)
            print(f'N-Ads bond lenght: {np.average(NAdsBondValues):.3f} Angs')
        for Ads_idx in mol.symbols.indices()[Ads]:
            Ads_position = mol[Ads_idx].position[2]
        print(f'Adsorbate position: {Ads_position:.3f}')
        print(f'Adsorbate distance from sheet: {abs(z_avg - Ads_position):.3f}')

    if show_mol: view(mol)

    # create a figure
    fig, axs = plt.subplots(2, 2)
    fig_full, axs_full = plt.subplots(2)

    # basic plot
    print(min_epot_idx, min_fmax_idx)
    idx_range = 5
    if min_epot_idx - idx_range < 0:
        epot_plot_idx = min_epot_idx + abs(min_epot_idx - idx_range)
    elif len(epot) - (min_epot_idx + idx_range) < 0:
        epot_plot_idx = abs(min_epot_idx - idx_range) + 1
    else:
        epot_plot_idx = min_epot_idx

    if min_fmax_idx - idx_range < 0:
        fmax_plot_idx = min_fmax_idx + abs(min_fmax_idx - idx_range)
    elif len(fmax) - (min_fmax_idx + idx_range) < 0:
        fmax_plot_idx = abs(min_fmax_idx - idx_range) + 1
    else:
        fmax_plot_idx = min_fmax_idx
    print(epot_plot_idx, fmax_plot_idx)

    if show_graph:
        # full plot
        axs_full[0].set_title('steps X epot (eV)')
        axs_full[0].plot(X, epot, color='blue')
        axs_full[1].set_title('steps X fmax (eV/A)')
        axs_full[1].plot(X, fmax, color='orange')

        # center around min energy
        axs[0, 0].set_title('steps X epot (eV)')
        axs[0, 0].plot(X[epot_plot_idx - idx_range:epot_plot_idx + idx_range],
                       epot[epot_plot_idx - idx_range:epot_plot_idx + idx_range],
                       'o-', color='blue')
        axs[1, 0].set_title('steps X fmax (eV/A)')
        axs[1, 0].plot(X[epot_plot_idx - idx_range:epot_plot_idx + idx_range],
                       fmax[epot_plot_idx - idx_range:epot_plot_idx + idx_range],
                       'o-', color='orange')

        # mark into the least energy point
        axs[0, 0].plot(X[min_epot_idx], epot[min_epot_idx], 'o', color='green')
        axs[1, 0].plot(X[min_epot_idx], fmax[min_epot_idx], 'o', color='green')
        axs_full[0].plot(X[min_epot_idx], epot[min_epot_idx], 'o', color='green')
        axs_full[1].plot(X[min_epot_idx], fmax[min_epot_idx], 'o', color='green')

        axs[0, 0].text(X[min_epot_idx], epot[min_epot_idx], str(epot[min_epot_idx]))
        axs[1, 0].text(X[min_epot_idx], fmax[min_epot_idx], str(fmax[min_epot_idx]))

        # center around min fmax
        axs[0, 1].set_title('steps X epot (eV)')
        axs[0, 1].plot(X[fmax_plot_idx - idx_range:fmax_plot_idx + idx_range],
                       epot[fmax_plot_idx - idx_range:fmax_plot_idx + idx_range],
                       'o-', color='blue')
        axs[1, 1].set_title('steps X fmax (eV/A)')
        axs[1, 1].plot(X[fmax_plot_idx - idx_range:fmax_plot_idx + idx_range],
                       fmax[fmax_plot_idx - idx_range:fmax_plot_idx + idx_range],
                       'o-', color='orange')

        # mark the least energy point
        axs[0, 1].plot(X[min_fmax_idx], epot[min_fmax_idx], 'o', color='red')
        axs[1, 1].plot(X[min_fmax_idx], fmax[min_fmax_idx], 'o', color='red')
        axs_full[0].plot(X[min_fmax_idx], epot[min_fmax_idx], 'o', color='red')
        axs_full[1].plot(X[min_fmax_idx], fmax[min_fmax_idx], 'o', color='red')

        axs[0, 1].text(X[min_fmax_idx], epot[min_fmax_idx], str(epot[min_fmax_idx]))
        axs[1, 1].text(X[min_fmax_idx], fmax[min_fmax_idx], str(fmax[min_fmax_idx]))

        plt.show()


