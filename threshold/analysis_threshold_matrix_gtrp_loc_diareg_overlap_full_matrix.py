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
    # import wide format data
    # df = pd.read_csv('gtrp_loc_diareg_overlap.txt', sep = '\t')

    # transform to long format
    # wide_to_long.expand(df, 'gtrp_loc_diareg_overlap_long.txt')

    df = pd.read_csv('gtrp_loc_diareg_overlap_long.txt', sep='\t')

    # remove diacritics
    df.transcription = [leven.remove_accents(trs) for trs in df.transcription]


    # generate character-numpy mappings and segment-segment cost matrix
    char_map, dec_dist, cost_mat = leven.init_cost_matrix_weighted('pmi_wieling_large.tsv')
    df.transcription  = leven.encode_column_trs_array(df.transcription, char_map)
    map_locations = {loc: idx for idx, loc in enumerate(df.location.unique())}

    # initialize location by location distance matrix
    result_mat = np.zeros((len(map_locations), len(map_locations)))

    for perm in permutations(map_locations, 2):
        loc1 = df[df.location == perm[0]]
        loc2 = df[df.location == perm[1]]

        distances = [leven.leven_compute_align(comp[0], comp[1], cost_mat)[0] for comp in zip(loc1.transcription, loc2.transcription) if comp[0][0] == -1 or comp[1][0] == -1]
        result_mat[map_locations[perm[0]], map_locations[perm[1]]] = mean_of_list(np.array(distances))

    np.savetxt('full_matrix_gtrp_loc_overlap_diareg.txt', result_mat, fmt='%.5f', delimiter='\t')
