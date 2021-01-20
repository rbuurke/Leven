from itertools import permutations

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


if __name__ == '__main__':
    # import RND data
    rnd_trs = pd.read_csv('rnd_long.txt', sep='\t')

    # generate character-numpy mappings and segment-segment cost matrix
    char_map, dec_dist, cost_mat = leven.init_cost_matrix_weighted('pmi_wieling.tsv')
    trs_enc_map = leven.generate_trs_map_df(rnd_trs.drop(columns=['location', 'word']), char_map)
    map_locations = {loc: idx for idx, loc in enumerate(rnd_trs.location.unique())}

    # initialize location by location distance matrix
    result_mat = np.zeros((len(map_locations), len(map_locations)))

    for perm in permutations(map_locations, 2):
        loc1 = rnd_trs[rnd_trs.location == perm[0]]
        loc2 = rnd_trs[rnd_trs.location == perm[1]]

        distances = [leven.leven_compute_align(trs_enc_map[comp[0]],
                                               trs_enc_map[comp[1]],
                                               cost_mat)[0] for comp in zip(loc1.transcription, loc2.transcription)
                     if comp[0] in trs_enc_map and comp[1] in trs_enc_map]

        result_mat[map_locations[perm[0]], map_locations[perm[1]]] = mean_of_list(np.array(distances))

    np.savetxt('full_matrix_rnd.txt', result_mat, fmt='%.5f', delimiter='\t')
