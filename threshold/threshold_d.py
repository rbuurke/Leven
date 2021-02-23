import os.path
import sys
import random
from itertools import permutations
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import pandas as pd
from numba import njit

sys.path.insert(0, '..')
import leven


@njit(cache=True)
def update_matrix(idx1, idx2, val, matrix):
    matrix[idx1, idx2] = val
    return matrix


def cycle_dict(computation_type, loc_ssize,
               word_ssize, transcriptions_dict,
               loc_mapping, dict_of_constants,
               list_of_variation, original_matrix):
    # print(loc_ssize, word_ssize)
    # initialize location by location distance matrix
    result_mat = np.zeros(original_matrix.shape)

    if computation_type == 'location':
        location_sample = dict_of_constants[loc_ssize]
        sample_of_words = sorted(
            random.sample(list_of_variation, word_ssize))

    if computation_type == 'word':
        location_sample = sorted(
            random.sample(list_of_variation, loc_ssize))
        sample_of_words = dict_of_constants[word_ssize]

    # compute distances between transcriptions of paired locations
    for perm in permutations(location_sample, 2):
        distances = []
        for word in sample_of_words:
            trs1 = transcriptions_dict[perm[0]][word]
            trs2 = transcriptions_dict[perm[1]][word]
            if not isinstance(trs1, float) \
                    and not isinstance(trs2, float):
                distances.append(leven.leven_dist(
                    trs1, trs2, cost_mat)[0])
        result_mat = update_matrix(loc_mapping[perm[0]],
                                   loc_mapping[perm[1]],
                                   np.mean(distances),
                                   result_mat)

    # select indices of the sampled locations and correlate this subset
    select_cols = sorted([map_locations[location]
                          for location in location_sample])

    vec_sample = result_mat[select_cols, :][:, select_cols]
    vec_orig = full_mat[select_cols, :][:, select_cols]

    # ensure that the diagonals are NaN and return correlation
    np.fill_diagonal(vec_orig, np.nan)
    np.fill_diagonal(vec_sample, np.nan)
    return loc_ssize, word_ssize, pd.Series(vec_orig.flatten()).corr(pd.Series(vec_sample.flatten()))


def meta_cycle(computation_type,
               dict_of_constants,
               list_of_variation):
    results = {cycle_dict(computation_type, loc_ssize,
                          word_ssize, df_dict,
                          map_locations, dict_of_constants,
                          list_of_variation, full_mat)
               for word_ssize in range(5, len(list_of_words), 5)
               for loc_ssize in dict_of_constants}
    print("Cycle done")
    return results


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
    char_map, dec_dist, cost_mat = leven.init_cost_matrix_weighted(
        '../pmi/pmi_wieling_large.tsv')
    df.transcription = leven.encode_column_trs_array_native(
        df.transcription, char_map)
    map_locations = {loc: idx for idx, loc in enumerate(df.location.unique())}

    # get all locations and words
    list_of_locs = list(df.location.unique())
    list_of_words = list(df.word.unique())

    random.seed(341995)

    ''' select random subsets: locations, variable words '''
    df_dict = {loc: {word: df.transcription[(df.word == word) & (df.location == loc)].iloc[0]
                     for word in df.word.unique()}
               for loc in df.location.unique()}
    #
    # loc_sample_df = {sample: random.sample(list_of_locs, sample)
    #                  for sample in range(5, len(list_of_locs), 5)}
    #
    # with ProcessPoolExecutor(max_workers=len(os.sched_getaffinity(0))) as executor:
    #     results_loc = [executor.submit(meta_cycle, 'location', loc_sample_df, list_of_words) for i in range(100)]
    #
    # output_df = pd.DataFrame([result for cycle in results_loc for result in cycle.result()],
    #                          columns=['loc_ssize', 'word_ssize', 'corr'])
    # output_df.sort_values(by=['loc_ssize', 'word_ssize']).to_csv('threshold_output_diareg_var_word.txt',
    #                                                              sep='\t',
    #                                                              index=False)

    ''' select random subsets: words, variable locations '''

    # df_dict = {loc: {word: df.transcription[(df.word == word) & (df.location == loc)].iloc[0]
    #                  for word in df.word.unique()}
    #            for loc in df.location.unique()}
    #
    # word_sample_df = {sample: random.sample(list_of_words, sample)
    #                   for sample in range(5, len(list_of_words), 5)}
    #
    # with ProcessPoolExecutor(max_workers=len(os.sched_getaffinity(0))) as executor:
    #     results_loc = [executor.submit(meta_cycle, 'word', word_sample_df, list_of_locs) for i in range(100)]
    #
    # output_df = pd.DataFrame([result for cycle in results_loc for result in cycle.result()],
    #                          columns=['loc_ssize', 'word_ssize', 'corr'])
    # output_df.sort_values(by=['loc_ssize', 'word_ssize']).to_csv('threshold_output_diareg_var_loc.txt',
    #                                                              sep='\t',
    #                                                              index=False)
