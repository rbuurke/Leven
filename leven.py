""" compute VC sensitive Levenshtein distance cost matrix
    retrace cost matrix to find VC sensitive alignment"""

import numpy as np
import pandas as pd
from numba import int64, float64
from numba import njit
from numba.core.types.containers import UniTuple
from numba.typed import List

list_type = UniTuple(int64, 2)

# set a limit on segment generation
alignments_lim = 50


def get_char_inv(df, drop_col=None):
    """ Get the character inventory of a DataFrame;
        Possible to remove a location column with the drop_col argument"""
    if drop_col:
        return {char for trs in df.drop(columns=drop_col).stack()
                for char in trs}
    else:
        return {char for trs in df.stack() for char in trs}


def remove_accents(input_str: str):
    """ Removes known diacritics from unicode type strings;
        Note that currently the '/' and spaces are also removed,
        although they are used to denote multiple transcriptions
        in Gabmap style data"""
    if isinstance(input_str, str):

        # if two transcriptions are given (in Gabmap format), take the first
        if '/' in input_str:
            input_str = input_str.split('/')[0]

        remove_chars_manual = {' ', '/', ':', '´', 'ʰ', 'ʲ', 'ʷ', 'ˀ',
                               'ˈ', 'ː', 'ˑ', '˞', 'ˠ', 'ˡ', 'ˢ', 'ˤ',
                               '˳', '˺', '̂', '̃', '̆', '̌', '̍', '̑',
                               '̘', '̙', '̚', '̜', '̝', '̞', '̟', '̠',
                               '̤', '̥', '̩', '̪', '̬', '̮', '̯', '̰',
                               '̱', '̹', '̺', '̻', '͂', '͆', '͉', '͎',
                               '͡', 'ᵊ', '̴', 'ʲ'}

        # manually remove chars
        final_form = ''.join(
            char for char in input_str if char not in remove_chars_manual)
        return final_form
    else:
        return input_str


def generate_trs_map(list_of_trs, char_map):
    """ Generates NumPy number encodings given a list """
    return {trs: np.array([char_map[char] for char in trs])
            for trs in list_of_trs if isinstance(trs, str)}


def encode_df(df, char_map):
    """ Generates NumPy number encodings given a dataframe """
    return {trs: np.array([char_map[char] for char in trs])
            for trs in df.stack() if isinstance(trs, str)}


def generate_char_map(list_of_chars):
    """ generate a character mapping from list of characters """
    enc_dict = {char: index for index, char in enumerate(list_of_chars)}
    dec_dict = dict(zip(enc_dict.values(), enc_dict.keys()))
    return enc_dict, dec_dict


def encode_column_trs_array(list_of_trs, char_map):
    """ encode a Pandas DataFrame column of transcriptions to NumPy arrays """
    encoded_list = []
    for trs in list_of_trs:
        if isinstance(trs, float):
            encoded_list.append(np.array([-1]))
        else:
            encoded_list.append(np.array([char_map[char] for char in trs]))
    return encoded_list


def encode_column_trs_array_native(list_of_trs, char_map):
    """ encode a Pandas DataFrame column of transcriptions to NumPy arrays """
    encoded_list = []
    for trs in list_of_trs:
        if isinstance(trs, float):
            encoded_list.append(np.nan)
        else:
            encoded_list.append(np.array([char_map[char] for char in trs]))
    return encoded_list


def encode_column_trs(list_of_trs, char_map):
    """ encode a Pandas DataFrame column of transcriptions to NumPy arrays """
    encoded_list = []
    for trs in list_of_trs:
        if isinstance(trs, float):
            encoded_list.append(np.nan)
        else:
            encoded_list.append(np.array([char_map[char] for char in trs]))
    return encoded_list


def init_cost_matrix_semivowels(list_of_chars, symbol_specification_dict):
    """ Generates a NumPy segment-segment matrix  **without** an input file
        Semivowels are allowed to match both vowels and consonants"""
    inv_size = len(list_of_chars)
    cost_matrix = np.zeros((inv_size, inv_size))

    for (i, j) in np.ndindex(inv_size, inv_size):
        # character inventory indices are used for filling in the cost matrix
        c1 = list_of_chars[i]
        c2 = list_of_chars[j]

        if c1 == c2:
            # match
            cost_matrix[i, j] = 0
        elif c1 == '-' or c2 == '-':
            # indel
            if c1 == '-':
                cost_matrix[-1, j] = 0.01
            if c2 == '-':
                cost_matrix[i, -1] = 0.01
        elif c1 in {'j', 'w'} or c2 in {'j', 'w'}:
            cost_matrix[i, j] = 0.01
        elif c1 != c2 and \
                symbol_specification_dict[c1] == symbol_specification_dict[c2]:
            # substitution
            cost_matrix[i, j] = 0.01
        elif c1 != c2 and \
                symbol_specification_dict[c1] != symbol_specification_dict[c2]:
            # V-C substitution
            cost_matrix[i, j] = 1
    return cost_matrix


def init_cost_matrix(list_of_chars, symbol_specification_dict):
    """ Generates a NumPy segment-segment matrix  **without** an input file """
    inv_size = len(list_of_chars)
    cost_matrix = np.zeros((inv_size, inv_size))

    for (i, j) in np.ndindex(inv_size, inv_size):
        # character inventory indices are used for filling in the cost matrix
        c1 = list_of_chars[i]
        c2 = list_of_chars[j]

        if c1 == c2:
            # match
            cost_matrix[i, j] = 0
        elif c1 == '-' or c2 == '-':
            # indel
            if c1 == '-':
                cost_matrix[-1, j] = 0.01
            if c2 == '-':
                cost_matrix[i, -1] = 0.01
        elif c1 != c2 and \
                symbol_specification_dict[c1] == symbol_specification_dict[c2]:
            # substitution
            cost_matrix[i, j] = 0.01
        elif c1 != c2 and \
                symbol_specification_dict[c1] != symbol_specification_dict[c2]:
            # V-C substitution
            cost_matrix[i, j] = 1
    return cost_matrix


def init_cost_matrix_weighted(pmi_file):
    """ Generates a NumPy segment-segment matrix  **with** an input file """

    with open(pmi_file, 'r') as file:
        char_list = file.readlines()[0].strip().split('\t')
        enc_dict, dec_dict = generate_char_map(char_list)

    dist_mat = pd.read_csv(pmi_file, sep='\t').drop(
        columns='Unnamed: 0').to_numpy()
    return enc_dict, dec_dict, dist_mat


@njit(cache=True)
def backtrack(a, array, w1_array, w2_array, alignment):
    """ backtrack the path from the pointer array without strings;
        faster for PMI matrix computation,
        although strings still need to be decoded for human legibility """
    # cut off alignment generation of total # of segments > set limit
    if len(alignment) > alignments_lim:
        return 0

    # paths = 0

    # for i in a[2:]:
    #     if not i == -1:
    #         paths += 1

    # if paths > 2:
    #     copy_alignment = [i for i in alignment]
    #     for i in copy_alignment:
    #         alignment.append(i)

    # if any of the targets are [0,0], stop tracing back
    if (a[2] == 0 and a[3] == 0) \
            or (a[4] == 0 and a[5] == 0) \
            or (a[6] == 0 and a[7] == 0):
        if a[2] == 0 and a[3] == 0:
            alignment.append((w1_array[a[2]], -1))
        if a[4] == 0 and a[5] == 0:
            alignment.append((-1, w2_array[a[5]]))
        if a[6] == 0 and a[7] == 0:
            alignment.append((w1_array[a[6]], w2_array[a[7]]))
        return 0

    for i, k in enumerate(a.reshape(4, 2)):
        # value is [-1, -1] or equal to current node, so continue iteration
        if np.max(k) == -1 \
                or (a[0] == k[0] and a[1] == k[1]):
            continue

        # append current node to alignment and do recursion
        else:
            if i == 1:
                # ins S1
                alignment.append((w1_array[a[2]], -1))

                # start again from connected nodes
                backtrack(array[(array[:, 0] == k[0]) &
                                (array[:, 1] == k[1])].flatten(),
                          array,
                          w1_array,
                          w2_array,
                          alignment)
            if i == 2:
                # ins S2
                alignment.append((-1, w2_array[a[5]]))

                # start again from connected nodes
                backtrack(array[(array[:, 0] == k[0]) &
                                (array[:, 1] == k[1])].flatten(),
                          array,
                          w1_array,
                          w2_array,
                          alignment)
            if i == 3:
                # substitute S1 and S2
                alignment.append((w1_array[a[6]], w2_array[a[7]]))

                # start again from connected nodes
                backtrack(array[(array[:, 0] == k[0]) &
                                (array[:, 1] == k[1])].flatten(),
                          array,
                          w1_array,
                          w2_array,
                          alignment)


@njit(cache=True)
def leven_compute_align(w1_idx, w2_idx, cost_matrix):
    """ compute the dynamic programming tableau;
        track 'pointers' from cell to cell"""
    row_nr = 0

    # init tableau
    tabl = np.zeros((w1_idx.size + 1, w2_idx.size + 1))

    # initiate the pointer matrix with lengths corresponding to input strings
    path_directions = np.full(
        (((w1_idx.size + 1) * (w2_idx.size + 1) - 1), 8), -1)

    # fill out first column and row + path directions
    for i, idx in enumerate(w1_idx):
        tabl[i + 1, 0] = tabl[i, 0] + cost_matrix[idx, -1]
        path_directions[row_nr] = np.array([i + 1, 0, i, 0, -1, -1, -1, -1])
        row_nr += 1

    for j, idx in enumerate(w2_idx):
        tabl[0, j + 1] = tabl[0, j] + cost_matrix[idx, -1]
        path_directions[row_nr] = np.array([0, j + 1, -1, -1, 0, j, -1, -1])
        row_nr += 1

    # begin iteration
    for (i, j) in np.ndindex(w1_idx.size + 1, w2_idx.size + 1):

        # ignore the first column and row
        if i != 0 and j != 0:
            path_directions[row_nr][0] = i
            path_directions[row_nr][1] = j

            # compute prior costs + added operation costs
            # going down: ins s1
            ins_s1 = tabl[i - 1, j] + cost_matrix[-1, w1_idx[i - 1]]

            # going right: ins s2
            ins_s2 = tabl[i, j - 1] + cost_matrix[-1, w2_idx[j - 1]]

            # going down-right
            sub = tabl[i - 1, j - 1] + cost_matrix[w1_idx[i - 1],
                                                   w2_idx[j - 1]]

            if (sub <= ins_s2) and (sub <= ins_s1):
                min_val = sub
                path_directions[row_nr][6] = i - 1
                path_directions[row_nr][7] = j - 1
            elif (ins_s2 <= sub) and (ins_s2 <= ins_s1):
                min_val = ins_s2
                path_directions[row_nr][4] = i
                path_directions[row_nr][5] = j - 1
            elif (ins_s1 <= sub) and (ins_s1 <= ins_s2):
                min_val = ins_s1
                path_directions[row_nr][2] = i - 1
                path_directions[row_nr][3] = j
            tabl[i, j] = min_val
            row_nr += 1
    alignment = List.empty_list(item_type=list_type)
    backtrack(path_directions[-1], path_directions, w1_idx, w2_idx, alignment)

    if len(alignment) > alignments_lim:
        return tabl[-1, -1], tabl, List.empty_list(item_type=list_type)
    else:
        return tabl[-1, -1], tabl, alignment

@njit(cache=True)
def leven_3_dim(w1_idx, w2_idx, cost_matrix):
    """ compute the dynamic programming tableau;
        track 'pointers' from cell to cell"""
    row_nr = 0

    # init tableau
    tabl = np.zeros((w1_idx.size + 1, w2_idx.size + 1))

    # initiate the pointer matrix with lengths corresponding to input strings
    path_directions = np.full(
        (((w1_idx.size + 1) * (w2_idx.size + 1) - 1), 8), -1)

    # fill out first column and row + path directions
    for i, idx in enumerate(w1_idx):
        tabl[i + 1, 0] = tabl[i, 0] + cost_matrix[idx, -1]
        path_directions[row_nr] = np.array([i + 1, 0, i, 0, -1, -1, -1, -1])
        row_nr += 1

    for j, idx in enumerate(w2_idx):
        tabl[0, j + 1] = tabl[0, j] + cost_matrix[idx, -1]
        path_directions[row_nr] = np.array([0, j + 1, -1, -1, 0, j, -1, -1])
        row_nr += 1

    # begin iteration
    for (i, j) in np.ndindex(w1_idx.size + 1, w2_idx.size + 1):

        # ignore the first column and row
        if i != 0 and j != 0:
            path_directions[row_nr][0] = i
            path_directions[row_nr][1] = j

            # compute prior costs + added operation costs
            # going down: ins s1
            ins_s1 = tabl[i - 1, j] + cost_matrix[-1, w1_idx[i - 1]]

            # going right: ins s2
            ins_s2 = tabl[i, j - 1] + cost_matrix[-1, w2_idx[j - 1]]

            # going down-right
            sub = tabl[i - 1, j - 1] + cost_matrix[w1_idx[i - 1],
                                                   w2_idx[j - 1]]

            if (sub <= ins_s2) and (sub <= ins_s1):
                min_val = sub
                path_directions[row_nr][6] = i - 1
                path_directions[row_nr][7] = j - 1
            elif (ins_s2 <= sub) and (ins_s2 <= ins_s1):
                min_val = ins_s2
                path_directions[row_nr][4] = i
                path_directions[row_nr][5] = j - 1
            elif (ins_s1 <= sub) and (ins_s1 <= ins_s2):
                min_val = ins_s1
                path_directions[row_nr][2] = i - 1
                path_directions[row_nr][3] = j
            tabl[i, j] = min_val
            row_nr += 1
    alignment = List.empty_list(item_type=list_type)
    backtrack(path_directions[-1], path_directions, w1_idx, w2_idx, alignment)

    if len(alignment) > alignments_lim:
        return tabl[-1, -1], tabl, List.empty_list(item_type=list_type)
    else:
        return tabl[-1, -1], tabl, alignment


@njit(cache=True)
def leven_dist(w1_idx, w2_idx, cost_matrix):
    """ compute the dynamic programming tableau"""
    if np.array_equal(w1_idx, w2_idx):
        # numba requires uniform return statements, hence these types
        return float64(0.0), np.zeros((1, 0))
    else:
        # init tableau
        tabl = np.zeros((w1_idx.size + 1, w2_idx.size + 1))

        # fill out first column and row + path directions
        for i, idx in enumerate(w1_idx):
            tabl[i + 1, 0] = tabl[i, 0] + cost_matrix[idx, -1]

        for j, idx in enumerate(w2_idx):
            tabl[0, j + 1] = tabl[0, j] + cost_matrix[idx, -1]

        # begin iteration
        for (i, j) in np.ndindex(w1_idx.size + 1, w2_idx.size + 1):

            # ignore the first column and row
            if i != 0 and j != 0:
                # going down: ins s1 / del s2
                ins_s1 = tabl[i - 1, j] + cost_matrix[-1, w1_idx[i - 1]]

                # going right: del s1 / ins s2
                ins_s2 = tabl[i, j - 1] + cost_matrix[-1, w2_idx[j - 1]]

                # going down-right
                sub = tabl[i - 1, j - 1] + cost_matrix[w1_idx[i - 1],
                                                       w2_idx[j - 1]]

                # compute prior costs + added operation costs
                costs = np.array([ins_s1, ins_s2, sub])

                # update tableau
                tabl[i, j] = np.min(costs)
        return tabl[-1, -1], tabl


if __name__ == '__main__':
    with open('inventory.txt', 'r', encoding='utf8') as inv_file:
        char_inv = {line.strip().split('\t')[0]: line.strip().split('\t')[1]
                    for line in inv_file.readlines()}
    chars_ = list(char_inv.keys())
    cost_mat = init_cost_matrix(chars_, char_inv)
    enc_map, dec_map = generate_char_map(chars_)

    dec_map[-1] = '-'

    ''' manual testing '''

    # str1 = 'kɑmərɑt'
    # str2 = 'vriəndkænɪsə'
    # str2 = 'vrind'
    # str1 = 'fbc'
    # str2 = 'ffbc'
    # str1 = 'ɑʊtos'
    # str2 = 'othos'
    # str1 = 'fabc'
    # str2 = 'abcd'
    # str1 = 'swart'
    # str2 = 'zwɑtɦ'
    str1 = 'abc'
    str2 = 'dbc'

    str1_ = np.array([enc_map[char] for char in str1])
    str2_ = np.array([enc_map[char] for char in str2])

    result = leven_compute_align(str1_, str2_, cost_mat)

    alignment_ = result[2]
    alignment_.reverse()

    decoded_alignment = [(dec_map[segment[0]], dec_map[segment[1]])
                         for segment in alignment_]

    # print(result[1])
    for segment in decoded_alignment:
        print(segment)
    # print(result[1])
    # print(result[3])

    ''' zwart testing '''
    # zwart_trs = pd.read_csv('trs_files/zwart_trs.tsv', sep='\t').transcription
    # for trs in zwart_trs:
    #     for comp in zwart_trs:
    #         if trs == comp or pd.isna(trs) or pd.isna(comp):
    #             continue
    #         trs_ = np.array([enc_map[char] for char in trs])
    #         comp_ = np.array([enc_map[char] for char in comp])

    #         result = leven_compute_align(trs_, comp_, cost_mat)
    #         print(trs, comp)
    #         alignment_ = result[2]
    #         alignment_.reverse()

    #         decoded_alignment = [(dec_map[segment[0]], dec_map[segment[1]])
    #                              for segment in alignment_]
    #         for segment in decoded_alignment:
    #             print(segment)
