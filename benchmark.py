import leven
from itertools import combinations
import numpy as np
import pandas as pd
from time import time

with open('inventory.txt', 'r', encoding='utf8') as inv_file:
    char_inv = {line.strip().split('\t')[0]: line.strip().split('\t')[1]
                for line in inv_file.readlines()}
chars_ = list(char_inv.keys())
cost_mat = leven.init_cost_matrix(chars_, char_inv)
enc_map, dec_map = leven.generate_char_map(chars_)
dec_map[-1] = '-'

# init the numba compilation
leven.leven_compute_align(np.array([enc_map[char] for char in 'abc']),
                          np.array([enc_map[char] for char in 'abc']), cost_mat)

# init the numba compilation
leven.leven_3_dim(np.array([enc_map[char] for char in 'abc']),
                  np.array([enc_map[char] for char in 'abc']),
                  np.array([enc_map[char] for char in 'abc']), cost_mat)

# import data
df = pd.read_csv('trs_files/GTRP.txt', sep='\t').iloc[:20, :]

# remove diacritics
for col in df.iloc[:, 1:]:
    df[col] = [leven.remove_accents(trs) for trs in df[col]]

# generate numpy mappings
enc_trs = {trs: np.array([enc_map[char] for char in trs])
           for col in df.iloc[:, 1:]
           for trs in df[col] if not pd.isna(trs)}


def cycle(column, data_frame):
    timings = []
    for i in range(10):
        start = time()
        count = 0
        data = []
        cache = {}
        for i, j in combinations(data_frame[column].dropna(), 2):
            if (i, j) in cache:
                data.append(cache[(i, j)])
            else:
                result = leven.leven_compute_align(enc_trs[i],
                                                   enc_trs[j],
                                                   cost_mat)
                cache[(i, j)] = result
            count += 1
        timings.append(time() - start)
    return count, sum(timings) / 10


def cycle_3d(column, data_frame):
    timings = []
    for i in range(10):
        start = time()
        count = 0

        data = []
        cache = {}
        for i, j, k in combinations(data_frame[column].dropna(), 3):
            if (i, j, k) in cache:
                data.append(cache[(i, j, k)])
            else:
                result = leven.leven_3_dim(enc_trs[i],
                                           enc_trs[j],
                                           enc_trs[k],
                                           cost_mat)
                cache[(i, j, k)] = result
            count += 1
        timings.append(time() - start)

    return count, sum(timings) / 10


for col in df.columns[1:5]:
    start = time()
    res = cycle(col, df)
    end = time()
    print(col, '\t', res[0], '\t', res[1])


for col in df.columns[1:5]:
    res = cycle_3d(col, df)
    print(col, '\t', res[0], '\t', res[1])
