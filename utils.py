def write_csv(out_path, file_name, **kwargs):
    import csv
    import numpy as np
    with open(f'{out_path}{file_name}.csv', 'w') as out_file:
        file_writer = csv.writer(out_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        cols_list = []
        for key, value in kwargs.items():
            # first I have to join key and value in a single list
            # after a little check, so I can work with arrays or list freely
            ar = np.array([])
            if type(value) == list:
                cols_list.append([key] + value)
            if type(value) == type(ar):
                value = value.tolist()
                cols_list.append([key] + value)
        rows = zip(*cols_list)
        for row in rows:
            file_writer.writerow(row)


def nearest_index(ref_list, number):
    import numpy as np
    arr = np.asarray(ref_list)
    i = (np.abs(arr - number)).argmin()
    return i


def path_dftb(path, dK, mol, verbose, get_dict):
    import numpy as np
    path2kpts = {'G': [0.0, 0.0, 0.0],
                 'M': [0.5, 0.0, 0.0],
                 'K': [2 / 3, 1 / 3, 0.0]}

    # I can't qute recall why this is here but is seems useful
    if get_dict:
        return path2kpts

    dftb_path = ''
    cell = mol.get_cell()
    reci = cell.reciprocal()
    k_f = path2kpts[path[0]]
    i = j = 0
    for point in path:
        k_i = k_f
        k_f = path2kpts[point]
        k_c = k_f[0] * reci[0] + k_f[1] * reci[1] + k_f[2] * reci[2]
        k_o = k_i[0] * reci[0] + k_i[1] * reci[1] + k_i[2] * reci[2]
        length = (k_c[0] - k_o[0]) ** 2 + (k_c[1] - k_o[1]) ** 2 + (k_c[2] - k_o[2]) ** 2
        length = np.sqrt(length)
        n = int(length / dK)
        if verbose >= 3:
            print(f'length from {path[j]} to {path[i]} is: {length}')
        if verbose >= 3:
            print(f'{n} steps of {dK} will add up to {n * dK}\n')
        if n == 0:
            n = 1
        dftb_path += f'{n}    {k_f[0]}    {k_f[1]}    {k_f[2]}\n'
        j = i
        i += 1
    return dftb_path

def read_detailed_dftb(path):
    with open(path) as detailed_file:
        for line in detailed_file:
            if 'Total energy:' in line:
                return line.split()[4]


def silly_method_to_get_plane_distance(molecule):
    # TODO: This is kinda stupid, I should do some proper math here
    import statistics
    positions = molecule.get_positions()
    z_ref = positions[0][2]
    z1_list = []
    z2_list = []
    for atom in positions:
        if abs(atom[2]-z_ref) > 1.0:
            z2_list.append(round(atom[2],3))
        else:
            z1_list.append(round(atom[2],3))
    #print(positions)
    #print(z1_list,z2_list)
    if z1_list and z2_list:
        pass
    elif not z1_list:
        z1_list = [round(positions[0][2],3)]
    elif not z2_list:
        z2_list = [round(positions[-1][2],3)]
    return abs(statistics.mean(z1_list) - statistics.mean(z2_list))

def read_dos_dftb(path, mol_name, return_eigen=False):
    ene = []
    dos = []
    i = 0
    print('read DOS')
    eigen_list = []
    with open(f'{path}{mol_name}.dos.dat') as file:
        for line in file:
            data = line.split()
            ene.insert(i, float(data[0]))
            dos.insert(i, float(data[1]))
            i += 1
    with open(f'{path}band.out') as eigen:
        for line in eigen:
            data = line.split()
            if len(data) == 3:
                eigen_list.append(float(data[1]))
    if return_eigen:
        return ene, dos, eigen_list
    else:
        return ene, dos


def read_orbital_energy_dftb(path, mol_name, orbitals):
    print(f'read {orbitals} orbitals')
    energies = []
    with open(f'{path}band.out') as inputFile:
        for line in inputFile:
            txt = line.split()
            if txt != [] and txt[0] != 'KPT':
                for orbital in orbitals:
                    if orbital == int(txt[0]):
                        energies.append(float(txt[1]))
    print(f'energies: {energies}')
    return energies


def read_fermi_levels_dftb(path, mol_name, verbosity=2):
    # POSSIBLE PROBLEM:
    # sometimes a Fermi Filling is required by some, but not all molecules in some set
    # if a zero gap material is present, then the input of a Fermi Filling certainly will ensure
    # some fraction of an occupation in valence band, making this search to find a "fake" half-filled band
    # this is problematic because I'd like to make Cody FLEXIBLE to run the same procedure for
    # various molecules, one possible avenue to solve this is do some rounding when searching for
    # occupation, certainly a 0.00001 occupation is irrelevant to many applications!
    print('read Fermi Level')
    lumo = [0, 0, 1000]
    homo = [0, 0, -1000]
    lumo_cur = [0, 0, 0]
    homo_cur = [0, 0, 0]
    with open(f'{path}band.out') as inputFile:
        for line in inputFile:
            txt = line.split()
            if txt != [] and txt[0] == 'KPT':  # true if we have header
                kpt = int(txt[1])  # useful for indirect gap?
            # true if NOT a header and not found lumo
            if txt != [] and txt[0] != 'KPT' and lumo_cur == [0, 0, 0]:
                band = int(txt[0])
                energy = float(txt[1])
                if float(txt[2]) < 1.0:  # found lowest unocupied KPOINT
                    lumo_cur = [kpt, band, energy]
                else:  # if occupation is not 0.0 then fill homo_cur
                    homo_cur = [kpt, band, energy]
                energy_prev = energy
            # true we have reached the end of a KPOINT info
            if not txt:  # empty list is a false!
                if lumo[2] >= lumo_cur[2]:
                    lumo = lumo_cur
                if homo[2] <= homo_cur[2]:
                    homo = homo_cur
                lumo_cur = [0, 0, 0]
                homo_cur = [0, 0, 0]
    fermi_e = round((lumo[2] + homo[2]) / 2, 6)
    gap = round(lumo[2] - homo[2], 6)
    if verbosity >=2:
        print(f'\n{"molName":<15} {"homo[kpt, Band, eV]":>20} {"lumo[kpt, Band, eV]":>20} {"gap":<3} {"fermi_e":<6}')
        homo_string = f'[{homo[0]},{homo[1]},{homo[2]}]'
        lumo_string = f'[{lumo[0]},{lumo[1]},{lumo[2]}]'
        print(f'{mol_name:<15} {homo_string:>20} {lumo_string:>20} {gap:<3} {fermi_e:<3}')
    return homo, lumo, gap, fermi_e


def output_fermi_levels_dftb(mol_name, homo, lumo, fermi_energy, gap):
    import os
    try:
        with open('FermiLevels.out', 'r') as inputFile:
            with open('FermiLevels.tmp', 'w+') as tmp:
                for line in inputFile:
                    # print(line)
                    tmp.write(line)
                next_line_list = [mol_name, str(homo), str(lumo), str(fermi_energy), str(gap)]
                next_line = '{: <50} {: >24} {: >24} {: >20} {: >20}'.format(*next_line_list)
                tmp.write(next_line + '\n')
        os.rename('FermiLevels.tmp', 'FermiLevels.out')
    except FileNotFoundError:
        with open('FermiLevels.out', 'w+') as inputFile:
            next_line_list = ['Molecule', 'HOMO', 'LUMO', 'FermiE', 'gap']
            next_line = '{: <50} {: >24} {: >24} {: >20} {: >20}'.format(*next_line_list)
            inputFile.write(next_line + '\n')
            next_line_list = [mol_name, str(homo), str(lumo), str(fermi_energy), str(gap)]
            next_line = '{: <50} {: >24} {: >24} {: >20} {: >20}'.format(*next_line_list)
            inputFile.write(next_line + '\n')
