import numpy as np


def convert2tab(fname):
    with open(fname) as file:
        char_num = int(file.readlines(1)[0].strip())

        output_matrix = np.zeros((char_num + 1, char_num + 1))
        input_rows = file.readlines()

        line_iter = 0

        for i in range(1, char_num + 1):
            for j in range(i):
                output_matrix[i, j] = float(input_rows[line_iter])
                line_iter += 1

    # make the matrix symmetrical
    output_matrix = output_matrix + output_matrix.T - \
        np.diag(np.diag(output_matrix))

    np.savetxt(fname.replace('dif', 'txt'),
               output_matrix, fmt='%.8f', delimiter='\t')
    return 0


convert2tab('tmp/sounddistances_pmi.dif')
convert2tab('tmp/sounddistances.dif')
