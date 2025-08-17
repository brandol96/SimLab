def read_array(line):
    line = line.split(':')[1]
    line = line.replace('[','')
    line = line.replace(']','')
    line = line.replace(' ','')
    line = line.split(',')
    line = [float(i) for i in line]
    return line

def read_z_IDX(line):
    print(line)
    line = line.split(':')[1]
    line = line.split(',')
    line = [float(i) for i in line]
    line[0] = int(line[0])
    return line[0]

def read_matrix_line(matrix_line):
    matrix_line = matrix_line.strip('\n')
    matrix_line = matrix_line.replace(' ','')
    matrix_line = matrix_line.replace('[','')
    matrix_line = matrix_line.replace(']','')
    matrix_line = matrix_line.split(',')
    if not matrix_line[-1]: matrix_line.pop()
    matrix_line = [float(i) for i in matrix_line]
    return matrix_line

def evaluate_grid(filename,verbose=False):
    import numpy as np
    import matplotlib
    read_flag = False
    energy_matrix = []
    energy_matrices = []
    z_idx_array = []
    with open(filename) as f:
        for line in f:
            if read_flag:
                if ']]' in line:
                    read_flag = False
                    energy_line = read_matrix_line(line)
                    energy_matrix.append(energy_line)
                    energy_matrices.append(energy_matrix)
                    energy_matrix = []
                else:
                    energy_line = read_matrix_line(line)
                    energy_matrix.append(energy_line)
            if 'X_ARRAY' in line:
                x_array = read_array(line)
            if 'Y_ARRAY' in line:
                y_array = read_array(line)
            if 'Z_ARRAY' in line:
                z_array = read_array(line)
            if 'Z_IDX' in line:
                z_idx_array.append(read_z_IDX(line))
                read_flag = True

    print('Read grid:')
    print('x: ',x_array,'\ny: ',y_array,'\nz: ', z_array)
    if verbose: print('z_idx: ',z_idx_array)
    min_en = 1e6
    min_z = min_x = min_y = 0
    for energy_matrix,z_idx in zip(energy_matrices,z_idx_array):
        energy_matrix = np.array(energy_matrix)
        if verbose: print(f'\n\nevaluating matrix: {z_idx} c = {z_array[z_idx]}')
        if verbose: print(energy_matrix, '\n')
        for lin in range(len(x_array)):
            for col in range(len(y_array)):
                x = x_array[lin]
                y = y_array[col]
                z = z_array[z_idx]
                E = energy_matrix[lin][col]
                if verbose: print(f'current point: x = {x}, y = {y}, c = {z}, E = {E}')
                if min_en > E:
                    min_en = E
                    min_x = x
                    min_y = y
                    min_z = z
    print('minimal energy configuration:')
    print(f'unit cell: ({min_x},{min_y},{min_z}) with energy: {min_en} eV')