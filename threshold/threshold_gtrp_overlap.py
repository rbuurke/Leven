import os.path
import sys
import glob
import random
from itertools import permutations
from multiprocessing import Pool
import numpy as np
import pandas as pd
import math
from numba import njit

sys.path.insert(0, '..')
import leven


@njit(cache=True)
def update_matrix(idx1, idx2, val, matrix):
    matrix[idx1, idx2] = val
    return matrix


def cycle_dict(sample_size_locations, sample_size_words,
               transcriptions_dict, mapping_of_locations,
               location_list, word_list,
               original_matrix, list_of_files):
    f_name = 'matrices/gtrp_diareg/gtrp_diareg_' \
             + str(sample_size_locations) \
             + '_word_' \
             + str(sample_size_words) \
             + '.txt'
    if f_name not in list_of_files:
        print(sample_size_locations, sample_size_words)
        correlations = []
        for bootstrap_sample in range(100):
            # initialize location by location distance matrix
            result_mat = np.zeros(original_matrix.shape)

            # draw random sample of locations and words
            location_sample = sorted(random.sample(location_list, sample_size_locations))
            sample_of_words = sorted(random.sample(word_list, sample_size_words))

            # compute distances between transcriptions of paired locations
            for perm in permutations(location_sample, 2):
                distances = []
                for word in sample_of_words:
                    trs1 = transcriptions_dict[perm[0]][word]
                    trs2 = transcriptions_dict[perm[1]][word]
                    if not isinstance(trs1, float) and not isinstance(trs2, float):
                        distances.append(leven.leven_dist(trs1, trs2, cost_mat)[0])
                result_mat = update_matrix(mapping_of_locations[perm[0]],
                                           mapping_of_locations[perm[1]],
                                           np.mean(distances),
                                           result_mat)

            # select the indices of the sampled locations and correlate this subset
            select_cols = sorted([map_locations[location] for location in location_sample])

            vec_sample = result_mat[select_cols, :][:, select_cols]
            vec_orig = full_mat[select_cols, :][:, select_cols]

            np.fill_diagonal(vec_orig, np.nan)
            np.fill_diagonal(vec_sample, np.nan)
            correlations.append(np.ma.corrcoef(np.ma.masked_invalid(vec_orig[np.tril_indices_from(vec_orig)]),
                                               np.ma.masked_invalid(vec_sample[np.tril_indices_from(
                                                   vec_sample)]))[0, -1])

        with open(f_name, 'w') as file:
            file.write('\t'.join([str(num) for num in correlations]))
    return 0


if __name__ == '__main__':
    ''' import matrix based on all transcriptions'''
    full_mat = np.loadtxt('gtrp_overlap_full_matrix.txt', delimiter='\t')

    ''' PMI weighted variant '''
    df = pd.read_csv('gtrp_loc_diareg_overlap_long.txt', sep='\t')

    # remove diacritics
    df.transcription = [leven.remove_accents(trs) for trs in df.transcription]

    # generate character-numpy mappings and segment-segment cost matrix
    char_map, dec_dist, cost_mat = leven.init_cost_matrix_weighted('../pmi/pmi_wieling_large.tsv')
    df.transcription = leven.encode_column_trs_array(df.transcription, char_map)
    map_locations = {loc: idx for idx, loc in enumerate(df.location.unique())}

    list_of_locs = list(df.location.unique())
    list_of_words = list(df.word.unique())

    file_list = set(glob.glob('matrices/gtrp_diareg/*.txt'))

    ''' select random subset '''
    random.seed(341995)

    df_dict = {loc: {word: df.transcription[(df.word == word) & (df.location == loc)].iloc[0]
                     for word in df.word.unique()}
               for loc in df.location.unique()}

    # the -1, because range() starts from 0 and sample must have a size of at least one less than population
    args = ((location_sample_size, word_sample_size, df_dict, map_locations, list_of_locs, list_of_words, full_mat,
             file_list) for
            location_sample_size in range(5, len(map_locations) - 1) for word_sample_size in range(5, 80 - 1))
    for i in args:
        cycle_dict(*i)

    # with Pool(processes=int(len(os.sched_getaffinity(0)) / 2) ) as p:
    #     data = p.starmap(cycle_dict, args)
    #     p.close()
    #     p.join()
