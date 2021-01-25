import os
import random
from concurrent.futures import ProcessPoolExecutor
from itertools import permutations

import numpy as np
import pandas as pd
from numba import njit
import sys

sys.path.insert(0, '..')
import leven

@njit(cache=True)
def mean_of_list(list_of_nr):
    if len(list_of_nr) > 0:
        return np.mean(list_of_nr)
    else:
        return np.nan


def compute_all_dist(trs_dict: dict, numpy_mapping: dict, permuted_list_locations):
    distances = {}
    for perm in permuted_list_locations:
        sum_list = np.array([leven.leven_compute_align(numpy_mapping[trs_dict[perm[0]][word]],
                                                       numpy_mapping[trs_dict[perm[1]][word]],
                                                       cost_mat)[0]
                             for word in trs_dict[perm[0]]
                             if trs_dict[perm[0]][word] in numpy_mapping and trs_dict[perm[1]][word] in numpy_mapping])
        distances.update({perm: mean_of_list(sum_list)})
    return distances


def cycle():
    dist_mat = np.zeros(shape=(len(df.location.unique()), len(df.location.unique())))

    perm_list = permutations(df.location.unique(), 2)
    selected_words = random.sample(list(df.word.unique()), sample_size)
    df_sel = df[df.word.isin(selected_words)].copy()

    df_dict = {loc: {trs_data.word: trs_data.transcription
                     for trs_data in df_sel[df_sel.location == loc].itertuples()}
               for loc in df_sel.location.unique()}

    results = compute_all_dist(df_dict, trs_enc_map, perm_list)

    for locations, mean in results.items():
        dist_mat[map_locations[locations[0]], map_locations[locations[1]]] = mean

    return dist_mat


if __name__ == '__main__':
    ''' import matrix based on all transcriptions'''
    full_mat = np.loadtxt('full_matrix_rnd.txt', delimiter='\t')

    ''' PMI weighted variant '''
    df = pd.read_csv('rnd_long.txt', sep='\t')

    char_map, dec_dist, cost_mat = leven.init_cost_matrix_weighted('pmi_wieling.tsv')
    trs_enc_map = leven.generate_trs_map_df(df.drop(columns=['location', 'word']), char_map)
    map_locations = {loc: idx for idx, loc in enumerate(df.location.unique())}

    ''' select random subset '''
    random.seed(341995)

    for sample_size in range(1, 101):
        print('Samples of size: ', sample_size)
        bootstrap_samples = 101

        with ProcessPoolExecutor(max_workers=len(os.sched_getaffinity(0))) as executor:
            data = {executor.submit(cycle) for bootstrap_sample in range(bootstrap_samples)}

        corr_list = [np.ma.corrcoef(full_mat.flatten(),
                                    np.ma.masked_invalid(result.result().flatten()))[0, -1]
                     for result in data]

        with open(str("aggregate_matrix_rnd_sample_" + str(sample_size) + ".txt"), 'w') as file:
            file.write('\t'.join([str(nr) for nr in corr_list]))
