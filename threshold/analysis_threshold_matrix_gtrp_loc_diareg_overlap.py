import os.path
import glob
import random
# from concurrent.futures import ProcessPoolExecutor
from itertools import permutations, repeat
from multiprocessing import Pool
import leven
import numpy as np
import pandas as pd
from numba import njit


@njit(cache=True)
def mean_of_list(list_of_nr):
    if len(list_of_nr) > 0:
        return np.mean(list_of_nr)
    else:
        return np.nan


def cycle(sample_size_locations, sample_size_words, transcriptions, mapping_of_locations, original_matrix, list_of_files):
    f_name = 'matrices/gtrp_diareg/gtrp_overlap_diareg_loc_' \
             + str(sample_size_locations) \
             + '_word_' \
             + str(sample_size_words) \
             + '.txt'

    if not f_name in list_of_files:
    # if os.path.exists(f_name):
        # print(sample_size_locations, sample_size_words)
        correlations = []
        for bootstrap_sample in range(100):
            # initialize location by location distance matrix
            result_mat = np.zeros(original_matrix.shape)

            # # draw random sample of locations and words
            location_sample = random.sample(list(transcriptions.location.unique()), sample_size_locations)
            sample_of_words = random.sample(list(transcriptions.word.unique()), sample_size_words)
            transcriptions_sel = transcriptions[
                (transcriptions.location.isin(location_sample)) & (transcriptions.word.isin(sample_of_words))]

            # location_sample = random.sample(list(transcriptions.location.unique()), sample_size_locations)
            # transcriptions_sel = transcriptions[transcriptions.location.isin(location_sample)]


            # sample_of_words = random.sample(list(transcriptions_sel.word.unique()), sample_size_words)
            # transcriptions_sel = transcriptions_sel[transcriptions_sel.word.isin(sample_of_words)]

            # compute distances between transcriptions of paired locations (using lazy evaluation to save memory)
            paired_transcriptions = [(perm[0], perm[1], np.mean(np.array(
                [leven.leven_dist(trs[0], trs[1], cost_mat)[0] for trs in
                 zip(transcriptions_sel[transcriptions_sel.location == perm[0]].transcription,
                     transcriptions_sel[transcriptions_sel.location == perm[1]].transcription) if
                 not trs[0][0] == -1 and not trs[1][0] == -1]))) for perm in permutations(location_sample, 2)]

            # update the distance matrix
            for result in paired_transcriptions:
                result_mat[mapping_of_locations[result[0]],
                           mapping_of_locations[result[1]]] = result[2]

            # select the indices of the sampled locations and correlate this subset
            select_cols = sorted(list({map_locations[location] for location in location_sample}))
            correlations.append(np.corrcoef(result_mat[select_cols, :][:, select_cols].flatten(),
                                            full_mat[select_cols, :][:, select_cols].flatten())[0, -1])
        with open(f_name, 'w') as file:
            file.write('\t'.join([str(num) for num in correlations]))
        return 0


if __name__ == '__main__':
    ''' import matrix based on all transcriptions'''
    full_mat = np.loadtxt('full_matrix_gtrp_loc_overlap_diareg.txt', delimiter='\t')

    ''' PMI weighted variant '''
    df = pd.read_csv('gtrp_loc_diareg_overlap_long.txt', sep='\t')

    # remove diacritics
    df.transcription = [leven.remove_accents(trs) for trs in df.transcription]

    # generate character-numpy mappings and segment-segment cost matrix
    char_map, dec_dist, cost_mat = leven.init_cost_matrix_weighted('../pmi/pmi_wieling_large.tsv')
    df.transcription = leven.encode_column_trs_array(df.transcription, char_map)
    map_locations = {loc: idx for idx, loc in enumerate(df.location.unique())}

    file_list = set(glob.glob('matrices/gtrp_diareg/*.txt'))

    ''' select random subset '''
    random.seed(341995)

    # the +1 for location sample size, because range() starts from 0
    args = ((location_sample_size + 1, word_sample_size, df, map_locations, full_mat, file_list) for location_sample_size in range(4, len(map_locations)) for word_sample_size in range(5, 101))

    with Pool(processes=len(os.sched_getaffinity(0))) as p:
        data = p.starmap(cycle, args)