def read_spec_ev(path, x_range):
    x = []
    y = []
    with open(f'{path}spec-ev.dat') as inp:
        for line in inp:
            data = line.split()
            if float(data[0]) >= x_range[0] and float(data[0]) <= x_range[1]:
                x.append(float(data[0]))
                y.append(float(data[1]))
    return x,y

def read_casida(out_path, cutoff_OscStr, cutoff_weight, energy_upper,
                peak_filter = False, output = False):
    # this will try to read an existing EXC.DAT inside the 'out_path' location
    # during this initial reading, the EXC.DAT will allways be filtered by
    # oscillator strenght and weight. It can be further filtered by providing a
    # peak_filter = [peak1,peak2,...] in units of eV to select a set of points in the
    # energy spectra to return the transitions arround the selected region
    # TODO: Implement something to make it possible to select the region around the peak_filter


    casida_list = []
    casida_full = []
    if output:
        print(f'energy (Ev),transition,weight,oscilator',
              file=open(f'{out_path}casida_filtered.out', 'w'))

    with open(f'{out_path}EXC.DAT') as casida:
        for line in casida:
            data = line.split()
            if '->' in data:
                energy = float(data[0])
                weight = float(data[5])
                oscilator = float(data[1])
                transition = [int(data[2]), int(data[4])]
                casida_full.append([energy, transition, weight, oscilator])
                if energy > energy_upper:
                    break
                else:
                    if weight >= cutoff_weight and oscilator >= cutoff_OscStr:
                        casida_list.append([energy, transition, weight, oscilator])
                        if output:
                            print(f'{energy},{transition},{weight},{oscilator}',
                                  file=open(f'{out_path}casida_filtered.out', 'a'))

    if peak_filter:
        casida_peak_filter = []
        for peak in peak_filter:
            print(f'target peak: {peak}')
            print(f'filter:{peak-0.05} | {peak+0.05}')
            for energy, transition, weight, oscilator in casida_list:
                if peak-0.05 <= energy <= peak + 0.05:
                    print(energy, transition, weight, oscilator)
                    casida_peak_filter.append([peak,energy, transition, weight, oscilator])

    if output:
        with open(f'{out_path}casida_peaks.out','w+') as out_file:
            for peak in peak_filter:
                i = 0
                for energy, transition, weight, oscilator in casida_full:
                    if energy > peak:
                        out_file.write(f'\npeak center: {peak:.3f}\n')
                        out_file.write('energy,transition,weight,oscilator\n')
                        aux = casida_full[i-3:i+3]
                        for entry in aux:
                            out_file.write(f'{entry[0]:.3f},{entry[1]},{entry[2]:.3f},{entry[3]:.8f}\n')
                        break
                    i += 1

    if peak_filter:
        return casida_peak_filter
    else:
        return casida_list