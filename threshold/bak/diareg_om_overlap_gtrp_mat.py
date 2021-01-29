from multiprocessing import Pool
import numpy as np
import pandas as pd
import os
import random
from itertools import permutations, combinations
import sys
sys.path.insert(0, '../..')
import leven
import wide_to_long

def cycle(loc1, loc2):
    return loc1, loc2, np.mean(np.array([leven.leven_dist(trs[0], trs[1], cost_mat)[0] for trs in
                                          zip(df[df.location == loc1].transcription,
                                              df[df.location == loc2].transcription) if
                                          not trs[0][0] == -1 and not trs[1][0] == -1]))


if __name__ == '__main__':
    # # import wide format data
    # df = pd.read_csv('diareg_om_gtrp.tsv', sep = '\t')

    # # transform to long format
    # wide_to_long.expand(df, 'diareg_om_gtrp_long.txt')

    df = pd.read_csv('../diareg_om_gtrp_long.txt', sep='\t')

    # remove diacritics
    df.transcription = [leven.remove_accents(trs) for trs in df.transcription]

    # ad hoc char replacement
    replace_dict = {'ɐ': 'ə', 'ʂ': 'ʃ', 'ɻ': 'ɹ'}
    df.transcription = df.transcription.replace(replace_dict, regex=True)

    # generate character-numpy mappings and segment-segment cost matrix
    char_map, dec_dist, cost_mat = leven.init_cost_matrix_weighted('../pmi/pmi_wieling_large.tsv')
    df.transcription = leven.encode_column_trs_array(df.transcription, char_map)
    map_locations = {loc: idx for idx, loc in enumerate(df.location.unique())}

    # initialize location by location distance matrix
    result_mat = np.zeros((len(map_locations), len(map_locations)))

    with Pool(processes=len(os.sched_getaffinity(0))) as p:
        data = p.starmap(cycle, [perm for perm in permutations(map_locations, 2)])
    
    for result in data:
        result_mat[map_locations[result[0]], map_locations[result[1]]] = result[2]

    np.savetxt('diareg_om_gtrp_full_matrix.txt', result_mat, fmt='%.5f', delimiter='\t')

    