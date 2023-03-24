from scipy import linalg, constants
import matplotlib.pyplot as plt
import numpy as np


def read_spec_ev(path, x_range):
    x = []
    y = []
    with open(f'{path}spec-ev.dat') as inp:
        for line in inp:
            data = line.split()
            if float(data[0]) >= x_range[0] and float(data[0]) <= x_range[1]:
                x.append(float(data[0]))
                y.append(float(data[1]))
    return x, y


def run_kick(method, out_path, mol_name, interactive_plot, directions, laser, fourier_damp, field_strength):
    # a few parameters
    zoom = 6
    title_font = 20
    label_font = 16
    text_font = 14
    x_range = [0, 30]

    print('\n\nstart plot\n\n')
    # setup figure
    fig = plt.figure()  # start a figure
    fig.suptitle(mol_name.replace("_", " "), fontsize=title_font)
    total = np.array([])
    time = np.array([])
    mu_matrix_is_empty = True

    d = 0.
    if 'X' in directions:
        print('\nLoading X direction...')
        mu = np.loadtxt(f'{out_path}mux.dat')  # response to excitation in x direction
        if total.size == 0:
            total = (mu[:, 1] - mu[0, 1])
        else:
            total += (mu[:, 1] - mu[0, 1])
        if time.size == 0:
            time = mu[:, 0]
        if mu_matrix_is_empty:
            mu_matrix = np.zeros((3, 3, mu.shape[0]))
            mu_matrix_is_empty = False
        mu_matrix[0, 0] = mu[:, 1]
        mu_matrix[0, 1] = mu[:, 2]
        mu_matrix[0, 2] = mu[:, 3]
        d += 1
    if 'Y' in directions:
        print('\nLoading Y direction...')
        mu = np.loadtxt(f'{out_path}muy.dat')  # response to excitation in y direction
        if total.size == 0:
            total = (mu[:, 2] - mu[0, 2])
        else:
            total += (mu[:, 2] - mu[0, 2])
        if time.size == 0:
            time = mu[:, 0]
        if mu_matrix_is_empty:
            mu_matrix = np.zeros((3, 3, mu.shape[0]))
            mu_matrix_is_empty = False
        mu_matrix[1, 0] = mu[:, 1]
        mu_matrix[1, 1] = mu[:, 2]
        mu_matrix[1, 2] = mu[:, 3]
        d += 1
    if 'Z' in directions:
        print('\nLoading Z direction...')
        mu = np.loadtxt(f'{out_path}muz.dat')  # response to excitation in z direction
        if total.size == 0:
            total = (mu[:, 3] - mu[0, 3])
        else:
            total += (mu[:, 3] - mu[0, 3])
        if time.size == 0:
            time = mu[:, 0]
        if mu_matrix_is_empty:
            mu_matrix = np.zeros((3, 3, mu.shape[0]))
            mu_matrix_is_empty = False
        mu_matrix[2, 0] = mu[:, 1]
        mu_matrix[2, 1] = mu[:, 2]
        mu_matrix[2, 2] = mu[:, 3]
        d += 1
    #some constants
    hplanck = constants.physical_constants['Planck constant in eV s'][0] * 1.0E15
    cspeednm = constants.speed_of_light * 1.0e9 / 1.0e15

    # full spectra calculation
    average = total / d
    damp = np.exp(-time / fourier_damp)
    field = field_strength

    spec = np.fft.rfft(damp * average, 10 * mu.shape[0])
    energsev = np.fft.rfftfreq(10 * mu.shape[0], mu[1, 0] - mu[0, 0]) * hplanck
    frec = np.fft.rfftfreq(10 * mu.shape[0], (mu[1, 0] - mu[0, 0]) * 1.0E-15)
    absorption = -2.0 * energsev * spec.imag / np.pi / field
    energsnm = constants.nu2lambda(frec[1:]) * 1.0E9

    emin = 0.5
    emax = 30.
    wvlmin = hplanck * cspeednm / emax
    wvlmax = hplanck * cspeednm / emin

    np.savetxt(f'{out_path}spec-ev.dat', np.column_stack((energsev[(energsev > emin) & (energsev < emax)], \
                                                          absorption[(energsev > emin) & (energsev < emax)])))
    np.savetxt(f'{out_path}spec-nm.dat', np.column_stack((energsnm[(energsnm > wvlmin) & (energsnm < wvlmax)], \
                                                          absorption[1:][(energsnm > wvlmin) & (energsnm < wvlmax)])))
    del energsev

    # read and plot the spectrum
    X, Y = read_spec_ev(out_path, x_range)
    plt.plot(X, Y, 'o-', markersize=2)

    # we want to evaluate where are the peaks to save their absorption and
    # calculate their maximum pol. direction
    N = len(Y)
    X_peaks = []
    Y_peaks = []
    for n in range(1, N - 1):
        if Y[n] > Y[n - 1] and Y[n] > Y[n + 1]:
            X_peaks.append(X[n])
            Y_peaks.append(Y[n])

    specs = np.zeros((3, 3, time.shape[0] * 5 + 1))
    alfa = np.zeros((3, 3))
    energsev = np.fft.rfftfreq(10 * time.shape[0], time[1] - time[0]) * hplanck
    X_pol = []
    Y_pol = []
    Z_pol = []
    for ene, abs in zip(X_peaks, Y_peaks):
        plt.plot([ene, ene], [0, abs], color='red')
        idx = (np.abs(energsev - ene)).argmin()
        for i in range(3):
            for j in range(3):
                specs[i, j] = - np.fft.rfft(damp * (mu_matrix[i, j, :] - mu_matrix[i, j, 0]),
                                            10 * time.shape[0]).imag * 2.0 * energsev / np.pi
                alfa[i, j] = specs[i, j, idx]
        w, v = linalg.eig(alfa)
        maxidx = np.argmax(w)
        pol_vector = (v[:,maxidx])
        X_pol.append(pol_vector[0])
        Y_pol.append(pol_vector[1])
        Z_pol.append(pol_vector[2])
        # print('PolarizationDirection = {:.8f} {:.8f} {:.8f}'.format(*v[:, maxidx]))

    # output csv using internal SimLab routine
    from SimLab.utils import write_csv
    file_name = f'{mol_name}_peaks_ev'
    write_csv(out_path, file_name, Energy=X_peaks, Absorption=Y_peaks,Pol_X=X_pol,Pol_Y=Y_pol,Pol_Z=Z_pol)

    plt.grid(True)
    fig.savefig(f'{out_path}{method}_{mol_name}_Optical.png')
    if interactive_plot:
        plt.show()
    else:
        plt.clf()


def run_laser(method, out_path, mol_name, interactive_plot, directions, laser, fourier_damp, field_strength):
    # a few parameters
    zoom = 6
    title_font = 20
    label_font = 16
    text_font = 14
    x_range = [0, 30]

    fig = plt.figure()  # start a figure
    fig.suptitle(mol_name.replace("_", " "), fontsize=title_font)
    mu = np.loadtxt(f'{out_path}mu.dat')
    T = []
    L = []
    for V in mu:
        T.append(round(V[0], 6))
        L.append(np.sqrt(V[1]**2 + V[2]**2 + V[3]**2))

    plt.plot(T, L, '-o', markersize=1.5)
    plt.grid(True)
    fig.savefig(f'{out_path}{method}_{mol_name}_laser.png')
    if interactive_plot:
        plt.show()
    else:
        plt.clf()