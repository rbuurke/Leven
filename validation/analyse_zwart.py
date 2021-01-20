import sys

import numpy as np
import pandas as pd

import leven


def compute_col(trs_list, trs_enc_map, cost_matrix):
    computations_col = []
    for target_trs in trs_list:
        for compare_trs in trs_list:
            if isinstance(target_trs, str) and isinstance(compare_trs, str):
                result = leven.leven_compute_align(trs_enc_map[target_trs],
                                                   trs_enc_map[compare_trs], cost_matrix)
                computations_col.append(result[0])
            else:
                computations_col.append(np.nan)
    return computations_col


np.set_printoptions(threshold=sys.maxsize)
df = pd.read_csv('zwart_trs.tsv', sep='\t')

''' PMI weighted variant '''
char_map, dec_dict, cost_mat = leven.init_cost_matrix_weighted('pmi_wieling.tsv')

trs_enc_map = leven.generate_trs_map(df.transcription, char_map)

results = compute_col(df.transcription, trs_enc_map, cost_mat)
# print(results)

output_mat = np.array(results).reshape((len(df.transcription), len(df.transcription)))

s1 = np.array([char_map[char] for char in 'zwat'])
s2 = np.array([char_map[char] for char in 'swɑt'])
dist = leven.leven_compute_align(s1, s2, cost_mat)
print(dist[1])

# print([cost_mat[char_map[char], -1] for char in 'zwart'])

np.savetxt('zwart_output_pmi.txt', output_mat, delimiter='\t', fmt='%.7f')
