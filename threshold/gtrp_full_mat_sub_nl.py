from multiprocessing import Pool
import numpy as np
import pandas as pd
import os
import random
from itertools import permutations, combinations
import sys

sys.path.insert(0, '..')
import leven


def cycle(loc1, loc2):
    return loc1, loc2, np.mean(np.array([leven.leven_dist(trs[0], trs[1], cost_mat)[0] for trs in
                                         zip(df[df.location == loc1].transcription,
                                             df[df.location == loc2].transcription) if
                                         not trs[0][0] == -1 and not trs[1][0] == -1]))


if __name__ == '__main__':
    # import wide format data
    # df = pd.read_csv('gtrp_loc_diareg_overlap.txt', sep = '\t')

    # transform to long format
    # wide_to_long.expand(df, 'gtrp_loc_diareg_overlap_long.txt')

    df = pd.read_csv('gtrp_loc_diareg_overlap_long.txt', sep='\t')

    # subset only Frisian, Groningen, and Drenthe locations
    with open('gtrp_provinces.txt', 'r') as file:
        map_locs = {line.strip().split('\t')[0].strip():
                    line.strip().split('\t')[1].strip()
                    for line in file.readlines()}

    filter_locs = [loc for loc in df.location.unique()
                   if map_locs[loc] in ['Groningen', 'Drenthe', 'Fryslân']]

    df = df[df.location.isin(filter_locs)]

    """ remove diacritics """
    df.transcription = [leven.remove_accents(trs) for trs in df.transcription]

    """ generate character-numpy mappings and segment-segment cost matrix """
    char_map, dec_dist, cost_mat = leven.init_cost_matrix_weighted('../pmi/pmi_wieling_large.tsv')
    df.transcription = leven.encode_column_trs_array(df.transcription, char_map)
    map_locations = {loc: idx for idx, loc in enumerate(df.location.unique())}

    """ initialize location by location distance matrix """
    result_mat = np.zeros((len(map_locations), len(map_locations)))

    args = ((loc1, loc2) for loc1 in map_locations for loc2 in map_locations)

    with Pool(processes=len(os.sched_getaffinity(0))) as p:
        data = p.starmap(cycle, [perm for perm in permutations(map_locations, 2)])

    for result in data:
        result_mat[map_locations[result[0]], map_locations[result[1]]] = result[2]

    np.savetxt('gtrp_overlap_full_matrix.txt', result_mat, fmt='%.5f', delimiter='\t')
