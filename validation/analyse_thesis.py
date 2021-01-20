import sys
import time
from itertools import permutations

import numpy as np
import pandas as pd

import leven


def compute_dist_col_dist_ptr_backtrack(list_of_trs, trs_enc_map, cache={}):
    # start iterating; 193 x 193 locations
    for trs, alt in permutations(list_of_trs, 2):
        if trs in trs_enc_map and alt in trs_enc_map:
            if (trs, alt) not in cache:
                result = leven.leven_compute_align(trs_enc_map[trs], trs_enc_map[alt], cost_mat)
                # leven.backtrack_str(result[2][-1], result[2], trs, alt)
                leven.backtrack(result[2][-1], result[2])
                cache[(trs, alt)] = result[0]
            else:
                cache[(trs, alt)]
        else:
            continue
    return 0


np.set_printoptions(threshold=sys.maxsize)
df = pd.read_csv('thesis_leven_test.tsv', sep='\t')

# remove diacritics
for col in df.drop(columns=['location']):
    df[col] = [leven.remove_accents(trs) if isinstance(trs, str) else trs for trs in df[col]]

''' PMI weighted variant '''
char_map, dec_dist, cost_mat = leven.init_cost_matrix_weighted('pmi_wieling.tsv')
trs_enc_map = leven.generate_trs_map_df(df.drop(columns=['location']), char_map)

# print(trs_enc_map)

''' tracking average times '''
timings = []

start_total = time.time()

for col in df.drop(columns=['location']):
    start = time.time()
    compute_dist_col_dist_ptr_backtrack(df[col], trs_enc_map)
    end = time.time()
    timings.append(end - start)
    print(col, 'finished')

end_total = time.time()
print('Average time per GTRP column: ', sum(timings) / len(timings), 's')
print('Total time: ', str(end_total - start_total), 's')
pd.Series(timings).to_csv('timings.txt', sep='\t', index=False)
