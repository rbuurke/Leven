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
    f_name = 'matrices/diareg_om_gtrp/diareg_om_gtrp_' \
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
    full_mat = np.loadtxt('diareg_om_gtrp_full_matrix.txt', delimiter='\t')

    ''' PMI weighted variant '''
    df = pd.read_csv('diareg_om_gtrp_long.txt', sep='\t')

    # remove diacritics
    df.transcription = [leven.remove_accents(trs) for trs in df.transcription]

    # ad hoc char replacement
    replace_dict = {'ɐ': 'ə', 'ʂ': 'ʃ', 'ɻ': 'ɹ'}
    df.transcription = df.transcription.replace(replace_dict, regex=True)

    # generate character-numpy mappings and segment-segment cost matrix
    char_map, dec_dist, cost_mat = leven.init_cost_matrix_weighted('../pmi/pmi_wieling_large.tsv')
    df.transcription = leven.encode_column_trs_array_native(df.transcription, char_map)
    map_locations = {loc: idx for idx, loc in enumerate(df.location.unique())}

    list_of_locs = list(df.location.unique())
    list_of_words = list(df.word.unique())

    file_list = set(glob.glob('matrices/diareg_om_gtrp/*.txt'))

    random.seed(341995)

    ''' select random subset '''

    df_dict = {loc: {word: df.transcription[(df.word == word) & (df.location == loc)].iloc[0]
                     for word in df.word.unique()}
               for loc in df.location.unique()}

    # the -1, because range() starts from 0 and sample must have a size of at least one less than population
    args = ((location_sample_size, word_sample_size, df_dict, map_locations, list_of_locs, list_of_words, full_mat,
             file_list) for
            location_sample_size in range(5, len(map_locations) - 1) for word_sample_size in range(5,
            len(list_of_words) - 1))
            # location_sample_size in range(5, len(map_locations) - 1, 5) for word_sample_size in range(5,
            # len(list_of_words) - 1, 5))

    with Pool(processes=len(os.sched_getaffinity(0))) as p:
        data = p.starmap(cycle_dict, args)
        p.close()
        p.join()

    # for i in args:
    #     cycle_dict(*i)

    ''' Testing correlations, recomputing full matrix with sample of words:
        20 words, adding one location at a time in dataframe order '''
    # df_dict = {loc: {word: df.transcription[(df.word == word) & (df.location == loc)].iloc[0]
    #                  for word in df.word.unique()}
    #            for loc in df.location.unique()}
    #
    # word_sample = random.sample(list_of_words, 20)
    #
    # # initiate list of correlations
    # correlations = []
    #
    # # initiate list of indices
    # locations = []
    #
    # # full matrix based on 20 words and all locations
    # full_mat_loc_sample = np.zeros(full_mat.shape)
    #
    # # compute distances between transcriptions of paired locations
    # for perm in permutations(list_of_locs, 2):
    #     distances = []
    #     for word in word_sample:
    #         trs1 = df_dict[perm[0]][word]
    #         trs2 = df_dict[perm[1]][word]
    #         if not isinstance(trs1, float) and not isinstance(trs2, float):
    #             distances.append(leven.leven_dist(trs1, trs2, cost_mat)[0])
    #     result_mat = update_matrix(map_locations[perm[0]],
    #                                map_locations[perm[1]],
    #                                np.mean(distances),
    #                                full_mat_loc_sample)
    # # print(full_mat_loc_sample)
    #
    # for loc_sample in list_of_locs:
    #     locations.append(loc_sample)
    #     if len(locations) > 1:
    #         # initialize location by location distance matrix
    #         result_mat = np.zeros(full_mat.shape)
    #
    #         # compute distances between transcriptions of paired locations
    #         for perm in permutations(locations, 2):
    #             distances = []
    #             for word in word_sample:
    #                 trs1 = df_dict[perm[0]][word]
    #                 trs2 = df_dict[perm[1]][word]
    #                 if not isinstance(trs1, float) and not isinstance(trs2, float):
    #                     distances.append(leven.leven_dist(trs1, trs2, cost_mat)[0])
    #             result_mat = update_matrix(map_locations[perm[0]],
    #                                        map_locations[perm[1]],
    #                                        np.mean(distances),
    #                                        result_mat)
    #
    #         # get indices for each location
    #         idx = [map_locations[loc] for loc in locations]
    #
    #         # select rows and columns from matrices
    #         vec_sample = result_mat[idx, :][:, idx]
    #         vec_orig = full_mat_loc_sample[idx, :][:, idx]
    #
    #
    #         correlations.append(np.corrcoef(vec_orig.flatten(), vec_sample.flatten())[0, -1])
    #     with open(str('matrices/diareg_om_gtrp/correlations_word_'
    #                   + str(len(word_sample))
    #                   + '_loc_'
    #                   + str(len(locations))
    #                   + '.txt'), 'w') as file:
    #         file.write('\n'.join([str(num) for num in correlations]))

    ''' Testing correlations, recomputing full matrix with sample of words:
        20 words, randomized sample order '''
    # df_dict = {loc: {word: df.transcription[(df.word == word) & (df.location == loc)].iloc[0]
    #                  for word in df.word.unique()}
    #            for loc in df.location.unique()}
    #
    # word_sample = random.sample(list_of_words, 20)
    #
    # # initiate list of correlations
    # correlations = []
    #
    # # full matrix based on 20 words and all locations
    # full_mat_loc_sample = np.zeros(full_mat.shape)
    #
    # # compute distances between transcriptions of paired locations
    # for perm in permutations(list_of_locs, 2):
    #     distances = []
    #     for word in word_sample:
    #         trs1 = df_dict[perm[0]][word]
    #         trs2 = df_dict[perm[1]][word]
    #         if not isinstance(trs1, float) and not isinstance(trs2, float):
    #             distances.append(leven.leven_dist(trs1, trs2, cost_mat)[0])
    #     result_mat = update_matrix(map_locations[perm[0]],
    #                                map_locations[perm[1]],
    #                                np.mean(distances),
    #                                full_mat_loc_sample)
    #
    # for sample_size in range(2, len(list_of_locs)):
    #     locations = random.sample(list_of_locs, sample_size)
    #
    #     # initialize location by location distance matrix
    #     result_mat = np.zeros(full_mat.shape)
    #
    #     # compute distances between transcriptions of paired locations
    #     for perm in permutations(locations, 2):
    #         distances = []
    #         for word in word_sample:
    #             trs1 = df_dict[perm[0]][word]
    #             trs2 = df_dict[perm[1]][word]
    #             if not isinstance(trs1, float) and not isinstance(trs2, float):
    #                 distances.append(leven.leven_dist(trs1, trs2, cost_mat)[0])
    #         result_mat = update_matrix(map_locations[perm[0]],
    #                                    map_locations[perm[1]],
    #                                    np.mean(distances),
    #                                    result_mat)
    #
    #     # get indices for each location
    #     idx = [map_locations[loc] for loc in locations]
    #
    #     # select rows and columns from matrices
    #     vec_sample = result_mat[idx, :][:, idx]
    #     vec_orig = full_mat_loc_sample[idx, :][:, idx]
    #
    #     correlations.append(np.corrcoef(vec_orig.flatten(), vec_sample.flatten())[0, -1])
    # with open(str('matrices/diareg_om_gtrp/correlations_word_'
    #               + str(len(word_sample))
    #               + '_loc_'
    #               + str(len(locations))
    #               + '.txt'), 'w') as file:
    #     file.write('\n'.join([str(num) for num in correlations]))
    # #
    # print(correlations)

    ''' Testing correlations, using original matrix based on *all* locations and words:
        20 words, linear order '''
    # df_dict = {loc: {word: df.transcription[(df.word == word) & (df.location == loc)].iloc[0]
    #                  for word in df.word.unique()}
    #            for loc in df.location.unique()}
    #
    # word_sample = sorted(random.sample(list_of_words, 80))
    # # word_sample = list_of_words
    #
    # # initiate list of correlations
    # correlations = []
    #
    # # initiate list of indices
    # locations = []
    #
    # for loc_sample in list_of_locs:
    #     locations.append(loc_sample)
    #     if 1 < len(locations):
    #         # initialize location by location distance matrix
    #         result_mat = np.zeros(full_mat.shape)
    #
    #         # compute distances between transcriptions of paired locations
    #         for perm in permutations(locations, 2):
    #             distances = []
    #             for word in word_sample:
    #                 trs1 = df_dict[perm[0]][word]
    #                 trs2 = df_dict[perm[1]][word]
    #                 if not isinstance(trs1, float) and not isinstance(trs2, float):
    #                     distances.append(leven.leven_dist(trs1, trs2, cost_mat)[0])
    #             result_mat = update_matrix(map_locations[perm[0]],
    #                                        map_locations[perm[1]],
    #                                        np.mean(distances),
    #                                        result_mat)
    #
    #         # get indices for each location
    #         idx = [map_locations[loc] for loc in locations]
    #
    #         # select rows and columns from matrices
    #         vec_sample = result_mat[idx, :][:, idx]
    #
    #         np.savetxt(str('matrices/diareg_om_gtrp/output_sample_mat_word_'
    #                        + str(len(word_sample))
    #                        + '_loc_'
    #                        + str(len(locations))
    #                        + '.txt'),
    #                    vec_sample,
    #                    fmt='%.5f')
    #
    #         vec_orig = full_mat[idx, :][:, idx]
    #
    #         np.savetxt(str('matrices/diareg_om_gtrp/output_original_mat_word_'
    #                        + str(len(word_sample))
    #                        + '_loc_'
    #                        + str(len(locations))
    #                        + '.txt'),
    #                    vec_orig,
    #                    fmt='%.5f')
    #         np.fill_diagonal(vec_orig, np.nan)
    #         np.fill_diagonal(vec_sample, np.nan)
    #         correlations.append(np.ma.corrcoef(np.ma.masked_invalid(vec_orig[np.tril_indices_from(vec_orig)]),
    #                                            np.ma.masked_invalid(vec_sample[np.tril_indices_from(vec_sample)]))[
    #                                 0, -1])
    #
    # with open(str('matrices/diareg_om_gtrp/correlations_word_linear_'
    #               + str(len(word_sample))
    #               + '_loc_'
    #               + str(len(list_of_locs))
    #               + '.txt'), 'w') as file:
    #     file.write('\n'.join([str(num) for num in correlations]))
    #
    # print(correlations)

    ''' Testing correlations, using original matrix based on *all* locations and words:
         20 words, randomized order '''
    # df_dict = {loc: {word: df.transcription[(df.word == word) & (df.location == loc)].iloc[0]
    #                  for word in df.word.unique()}
    #            for loc in df.location.unique()}
    #
    # word_sample = random.sample(list_of_words, 80)
    #
    # # initiate list of correlations
    # correlations = []
    #
    # # full matrix based on 20 words and all locations
    # full_mat_loc_sample = np.zeros(full_mat.shape)
    #
    # for sample_size in range(2, len(list_of_locs)):
    #     locations = random.sample(list_of_locs, sample_size)
    #
    #     # initialize location by location distance matrix
    #     result_mat = np.zeros(full_mat.shape)
    #
    #     # compute distances between transcriptions of paired locations
    #     for perm in permutations(locations, 2):
    #         distances = []
    #         for word in word_sample:
    #             trs1 = df_dict[perm[0]][word]
    #             trs2 = df_dict[perm[1]][word]
    #             if not isinstance(trs1, float) and not isinstance(trs2, float):
    #                 distances.append(leven.leven_dist(trs1, trs2, cost_mat)[0])
    #         result_mat = update_matrix(map_locations[perm[0]],
    #                                    map_locations[perm[1]],
    #                                    np.mean(distances),
    #                                    result_mat)
    #
    #     # get indices for each location
    #     idx = [map_locations[loc] for loc in locations]
    #
    #     # select rows and columns from matrices
    #     vec_sample = result_mat[idx, :][:, idx]
    #
    #     np.savetxt(str('matrices/diareg_om_gtrp/output_sample_mat_word_'
    #                    + str(len(word_sample))
    #                    + '_loc_'
    #                    + str(len(locations))
    #                    + '.txt'),
    #                vec_sample,
    #                fmt='%.5f')
    #
    #     vec_orig = full_mat[idx, :][:, idx]
    #
    #     np.savetxt(str('matrices/diareg_om_gtrp/output_original_mat_word_'
    #                    + str(len(word_sample))
    #                    + '_loc_'
    #                    + str(len(locations))
    #                    + '.txt'),
    #                vec_orig,
    #                fmt='%.5f')
    #
    #     np.fill_diagonal(vec_orig, np.nan)
    #     np.fill_diagonal(vec_sample, np.nan)
    #     correlations.append(np.ma.corrcoef(np.ma.masked_invalid(vec_orig[np.tril_indices_from(vec_orig)]),
    #                                        np.ma.masked_invalid(vec_sample[np.tril_indices_from(vec_sample)]))[0, -1])
    # with open(str('matrices/diareg_om_gtrp/correlations_word_'
    #               + str(len(word_sample))
    #               + '_loc_'
    #               + str(len(list_of_locs))
    #               + '.txt'), 'w') as file:
    #     file.write('\n'.join([str(num) for num in correlations]))
    #
    # print(correlations)
