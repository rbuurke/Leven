from itertools import permutations

import sys
import numpy as np
import pandas as pd
from numba import njit

sys.path.insert(0, '../..')
import leven


@njit(cache=True)
def mean_of_list(list_of_nr):
    if len(list_of_nr) > 0:
        return np.mean(list_of_nr)
    else:
        return np.nan


if __name__ == '__main__':
    # # import wide format data
    # df = pd.read_csv('transcription_files - DiaReg_om_stripped_diacritics_priority.tsv', sep = '\t')

    # # transform to long format
    # wide_to_long.expand(df, 'diareg_om_long.txt')


    # # import DiaReg OM data
    # df = pd.read_csv('diareg_om_long.txt', sep='\t')

    # ''' replace chars '''
    # replacement_dict = {'ɐ': 'ə', '̴': '', 'ʂ': 'ʃ', 'ɻ' : 'ɹ', 'ɬ': 'l'}

    # trs_fixed = []
    # for trs in df.transcription:
    #     if isinstance(trs, str):
    #         trs_fixed.append(''.join([replacement_dict[char] if char in replacement_dict else char for char in trs]))
    #     else:
    #         trs_fixed.append(trs)

    # df.transcription = trs_fixed

    # # export to new file and continue
    # df.to_csv('diareg_om_long_fixed.txt', sep='\t', index=False)

    df = pd.read_csv('diareg_om_long_fixed.txt', sep='\t')

    # generate character-numpy mappings and segment-segment cost matrix
    char_map, dec_dist, cost_mat = leven.init_cost_matrix_weighted('pmi_wieling.tsv')
    trs_enc_map = leven.generate_trs_map_df(df.drop(columns=['location', 'word']), char_map)

    map_locations = {loc: idx for idx, loc in enumerate(df.location.unique())}

    # initialize location by location distance matrix
    result_mat = np.zeros((len(map_locations), len(map_locations)))

    for perm in permutations(map_locations, 2):
        loc1 = df[df.location == perm[0]]
        loc2 = df[df.location == perm[1]]

        distances = [leven.leven_compute_align(trs_enc_map[comp[0]],
                                               trs_enc_map[comp[1]],
                                               cost_mat)[0] for comp in zip(loc1.transcription, loc2.transcription)
                     if comp[0] in trs_enc_map and comp[1] in trs_enc_map]

        result_mat[map_locations[perm[0]], map_locations[perm[1]]] = mean_of_list(np.array(distances))

    np.savetxt('full_matrix_diareg_om.txt', result_mat, fmt='%.5f', delimiter='\t')
