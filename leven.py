""" compute VC sensitive Levenshtein distance cost matrix
    retrace cost matrix to find VC sensitive alignment"""

import resource
from sys import maxsize, setrecursionlimit
import numpy as np
import pandas as pd
from numba import int64, float64
from numba import njit
from numba.core.types.containers import UniTuple
from numba.typed import List

list_type = UniTuple(int64, 2)
setrecursionlimit(10**6)


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
def member_(node, trace_list):
    is_visited = False
    for check_node in trace_list:
        if (node[0], node[1]) == check_node:
            is_visited = True
    return is_visited


@njit(cache=True)
def check_first(node, list_):
    if len(list_) == 0:
        return False
    else:
        if (node[0], node[1]) == list_[0]:
            return True
        else:
            return False


@njit(cache=True)
def add_alignment(current_alignment, aggregate):
    for node in current_alignment:
        aggregate.append(node)


@njit(cache=True)
def check_final(node, valid_vals):
    """ Check if current option is the final one"""
    if (node[0], node[1]) == (valid_vals[-2], valid_vals[-1]):
        return True
    else:
        return False


@njit(cache=True)
def check_depth(node, trace_list):
    if len(trace_list) > 0:
        if trace_list[0][0] < node[0] or trace_list[0][1] < node[1]:
            return True
        else:
            return False
    else:
        return False


@njit(cache=True)
def backtrack(a, array, w1_array, w2_array, trace_list, align_, alignment):
    """ backtrack the path from the pointer array without strings;
        faster for PMI matrix computation,
        although strings still need to be decoded for human legibility """
    # construct reshaped array of only source node coordinates
    iter_ = a.reshape(4, 2)[1:]

    # compute the number of options (note: value is doubled)
    options = np.sum(np.not_equal(iter_, -1))

    # view only valid coordinates for index checking (see below)
    valid_ = a[a != -1]
    # print(a)
    # print(trace_list)

    for i, k in enumerate(iter_):
        # value is [-1, -1], so continue iteration
        if np.max(k) == -1:
            continue
        # value is [0, 0], so we are at the beginning
        if np.max(k) == 0:
            if i == 0:
                align_.append((w1_array[a[2]], -1))
            if i == 1:
                align_.append((-1, w2_array[a[5]]))
            if i == 2:
                align_.append((w1_array[a[6]], w2_array[a[7]]))

            # align_.append((-1, -1))
            add_alignment(align_, alignment)
            if len(trace_list) == 0:
                return 0
            else:
                align_.clear()
                return backtrack(array[-1], array, w1_array, w2_array,
                                 trace_list, align_, alignment)
        else:
            if options > 2:
                final_ = check_final(k, valid_)

                if member_(k, trace_list):
                    if i == 0:
                        if check_first(k, trace_list):
                            trace_list.pop(0)
                        else:
                            align_.append((w1_array[a[2]], -1))
                            return backtrack(array[(array[:, 0] == k[0])
                                                   & (array[:, 1] == k[1])].
                                             flatten(),
                                             array, w1_array, w2_array,
                                             trace_list, align_,
                                             alignment)
                    if i == 1:
                        if check_first(k, trace_list):
                            trace_list.pop(0)
                        else:
                            align_.append((-1, w2_array[a[5]]))
                            return backtrack(array[(array[:, 0] == k[0])
                                                   & (array[:, 1] == k[1])].
                                             flatten(),
                                             array, w1_array, w2_array,
                                             trace_list, align_,
                                             alignment)
                    if i == 2:
                        if check_first(k, trace_list):
                            trace_list.pop(0)
                        else:
                            align_.append((w1_array[a[6]], w2_array[a[7]]))
                            return backtrack(array[(array[:, 0] == k[0])
                                                   & (array[:, 1] == k[1])].
                                             flatten(),
                                             array, w1_array, w2_array,
                                             trace_list, align_,
                                             alignment)
                else:
                    if check_depth(k, trace_list):
                        if final_:
                            if i == 0:
                                align_.append((w1_array[a[2]], -1))
                                return backtrack(array[(array[:, 0] == k[0])
                                                       & (array[:, 1] == k[1])].
                                                 flatten(),
                                                 array, w1_array, w2_array,
                                                 trace_list, align_,
                                                 alignment)
                            if i == 1:
                                align_.append((-1, w2_array[a[5]]))
                                return backtrack(array[(array[:, 0] == k[0])
                                                       & (array[:, 1] == k[1])].
                                                 flatten(),
                                                 array, w1_array, w2_array,
                                                 trace_list, align_,
                                                 alignment)
                            if i == 2:
                                align_.append((w1_array[a[6]], w2_array[a[7]]))
                                return backtrack(array[(array[:, 0] == k[0])
                                                       & (array[:, 1] == k[1])].
                                                 flatten(),
                                                 array, w1_array, w2_array,
                                                 trace_list, align_,
                                                 alignment)
                    else:
                        if final_:
                            if i == 0:
                                align_.append((w1_array[a[2]], -1))
                                return backtrack(array[(array[:, 0] == k[0])
                                                       & (array[:, 1] == k[1])].
                                                 flatten(),
                                                 array, w1_array, w2_array,
                                                 trace_list, align_,
                                                 alignment)
                            if i == 1:
                                align_.append((-1, w2_array[a[5]]))
                                return backtrack(array[(array[:, 0] == k[0])
                                                       & (array[:, 1] == k[1])].
                                                 flatten(),
                                                 array, w1_array, w2_array,
                                                 trace_list, align_,
                                                 alignment)
                            if i == 2:
                                align_.append((w1_array[a[6]], w2_array[a[7]]))
                                return backtrack(array[(array[:, 0] == k[0])
                                                       & (array[:, 1] == k[1])].
                                                 flatten(),
                                                 array, w1_array, w2_array,
                                                 trace_list, align_,
                                                 alignment)
                        else:
                            if i == 0:
                                trace_list.insert(0, (k[0], k[1]))
                                align_.append((w1_array[a[2]], -1))
                                return backtrack(array[(array[:, 0] == k[0])
                                                       & (array[:, 1] == k[1])].
                                                 flatten(),
                                                 array, w1_array, w2_array,
                                                 trace_list, align_,
                                                 alignment)
                            if i == 1:
                                trace_list.insert(0, (k[0], k[1]))
                                align_.append((-1, w2_array[a[5]]))
                                return backtrack(array[(array[:, 0] == k[0])
                                                       & (array[:, 1] == k[1])].
                                                 flatten(),
                                                 array, w1_array, w2_array,
                                                 trace_list, align_,
                                                 alignment)
                            if i == 2:
                                trace_list.insert(0, (k[0], k[1]))
                                align_.append((w1_array[a[6]], w2_array[a[7]]))
                                return backtrack(array[(array[:, 0] == k[0])
                                                       & (array[:, 1] == k[1])].
                                                 flatten(),
                                                 array, w1_array, w2_array,
                                                 trace_list, align_,
                                                 alignment)
            else:
                if i == 0:
                    align_.append((w1_array[a[2]], -1))
                    return backtrack(array[(array[:, 0] == k[0])
                                           & (array[:, 1] == k[1])].flatten(),
                                     array, w1_array, w2_array,
                                     trace_list, align_,
                                     alignment)

                if i == 1:
                    align_.append((-1, w2_array[a[5]]))
                    return backtrack(array[(array[:, 0] == k[0])
                                           & (array[:, 1] == k[1])].flatten(),
                                     array, w1_array, w2_array,
                                     trace_list, align_,
                                     alignment)
                if i == 2:
                    align_.append((w1_array[a[6]], w2_array[a[7]]))
                    return backtrack(array[(array[:, 0] == k[0])
                                           & (array[:, 1] == k[1])].flatten(),
                                     array, w1_array, w2_array,
                                     trace_list, align_,
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
            if (ins_s2 <= sub) and (ins_s2 <= ins_s1):
                min_val = ins_s2
                path_directions[row_nr][4] = i
                path_directions[row_nr][5] = j - 1
            if (ins_s1 <= sub) and (ins_s1 <= ins_s2):
                min_val = ins_s1
                path_directions[row_nr][2] = i - 1
                path_directions[row_nr][3] = j
            # tabl[i, j] = min(ins_s1, ins_s2, sub)
            tabl[i, j] = min_val
            row_nr += 1
    alignment = List.empty_list(item_type=list_type)
    current_ = List.empty_list(item_type=list_type)
    trace = List.empty_list(item_type=list_type)
    backtrack(path_directions[-1], path_directions, w1_idx, w2_idx, trace,
          current_, alignment)
    return tabl[-1, -1], tabl, alignment


@njit(cache=True)
def weight_(c1_idx, c2_idx, c3_idx, cost_matrix):
    """ compute the sum of all pairwise weights"""
    return cost_matrix[c1_idx, c2_idx] \
           + cost_matrix[c1_idx, c3_idx] \
           + cost_matrix[c2_idx, c3_idx]


@njit(cache=True)
def min_(o1, o2, o3, o4, o5, o6, o7):
    """ returns the minimum value from the given operation costs"""
    if o1 < o2 and o1 < o3 and o1 < o4 and o1 < o5 and o1 < o6 and o1 < o7:
        return o1
    if o2 < o3 and o2 < o4 and o2 < o5 and o2 < o6 and o2 < o7:
        return o2
    if o3 < o4 and o3 < o5 and o3 < o6 and o3 < o7:
        return o3
    if o4 < o5 and o4 < o6 and o4 < o7:
        return o4
    if o5 < o6 and o5 < o7:
        return o5
    if o6 < o7:
        return o6
    else:
        return o7


@njit(cache=True)
def leven_3_dim(w1_idx, w2_idx, w3_idx, cost_matrix):
    """ compute the dynamic programming tableau;
        track 'pointers' from cell to cell"""
    # row_nr = 0

    ''' tableau initialization
    note that the maximum length of the 3 strings is used here'''
    tabl = np.zeros((w3_idx.size + 1, w2_idx.size + 1, w1_idx.size + 1))

    # initiate the pointer matrix with lengths corresponding to input strings
    # path_directions = np.full(((tabl.size - 1), 24), -1)

    # begin iteration
    for (i, j, k), val in np.ndenumerate(tabl):
        # init the operation costs to a large number
        i1 = maxsize
        i2 = maxsize
        i3 = maxsize
        i12 = maxsize
        i13 = maxsize
        i23 = maxsize
        i123 = maxsize

        # determine costs for each operation
        if i > 0:
            i1 = tabl[i - 1, j, k] + weight_(-1, -1, w3_idx[i - 1],
                                             cost_matrix)
        if j > 0:
            i2 = tabl[i, j - 1, k] + weight_(-1, w2_idx[j - 1], -1,
                                             cost_matrix)
        if k > 0:
            i3 = tabl[i, j, k - 1] + weight_(w1_idx[k - 1], -1, -1,
                                             cost_matrix)
        if i > 0 and j > 0:
            i23 = tabl[i - 1, j - 1, k] + weight_(-1, w2_idx[j - 1],
                                                  w3_idx[i - 1],
                                                  cost_matrix)
        if i > 0 and k > 0:
            i13 = tabl[i - 1, j, k - 1] + weight_(w1_idx[k - 1], -1,
                                                  w3_idx[i - 1],
                                                  cost_matrix)
        if j > 0 and k > 0:
            i12 = tabl[i, j - 1, k - 1] + weight_(-1, w2_idx[j - 1],
                                                  w1_idx[k - 1],
                                                  cost_matrix)
        if i > 0 and j > 0 and k > 0:
            i123 = tabl[i - 1, j - 1, k - 1] + weight_(w1_idx[k - 1],
                                                       w2_idx[j - 1],
                                                       w1_idx[k - 1],
                                                       cost_matrix)
        min_val = min_(i1, i2, i3, i23, i13, i12, i123)
        tabl[i, j, k] = min_val
    return tabl


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
    # str1 = 'abc'
    # str2 = 'dbc'

    # str1_ = np.array([enc_map[char] for char in str1])
    # str2_ = np.array([enc_map[char] for char in str2])

    # result = leven_compute_align(str1_, str2_, cost_mat)
    #
    # alignment_ = result[2]
    # alignment_.reverse()
    #
    # decoded_alignment = [(dec_map[segment[0]], dec_map[segment[1]])
    #                      for segment in alignment_]

    # for segment in decoded_alignment:
    #     print(segment)

    ''' zwart testing '''
    # zwart_trs = pd.read_csv('trs_files/zwart_trs.tsv', sep='\t').transcription
    # for trs in zwart_trs:
    #     for comp in zwart_trs:
    #         print(trs, comp)
    #         if trs == comp or pd.isna(trs) or pd.isna(comp):
    #             continue
    #         trs_ = np.array([enc_map[char] for char in trs])
    #         comp_ = np.array([enc_map[char] for char in comp])
    #
    #         result = leven_compute_align(trs_, comp_, cost_mat)
    #         alignment_ = result[2]
    #         alignment_.reverse()
    #
    #         decoded_alignment = [(dec_map[segment[0]], dec_map[segment[1]])
    #                              for segment in alignment_]
    #         for segment in decoded_alignment:
    #             print(segment)

    ''' 3 dim testing '''
    # str1 = 'fabc'
    # str2 = 'aec'
    # str1 = 'bɪndən'
    # str1 = 'swart'
    # str2 = 'bɛində'
    # str1 = 'ʔæχt'
    # str2 = 'ouədə'
    str1 = 'ʔɒrde'
    str2 = 'ɛət'
    # str1 = 'kipən'
    # str2 = 'hɔundrɔundrhɔunrɔunr'
    str3 = 'abc'

    str1_ = np.array([enc_map[char] for char in str1])
    str2_ = np.array([enc_map[char] for char in str2])
    str3_ = np.array([enc_map[char] for char in str3])

    # result = leven_3_dim(str1_, str2_, str3_, cost_mat)
    # print(result)

    # for i in range(1 + 1):
    #     for j in range(2 + 1):
    #         for k in range(3 + 1):
    #             print(i, j, k)

    # test_tabl = np.arange(2*3*4)
    # print(test_tabl)

    # test_tabl = test_tabl.reshape((2, 3, 4))
    # print(test_tabl)
    # print((test_tabl[0, 0, 0]))

    # for i in np.ndenumerate(test_tabl):
    # print(i)

    result = leven_compute_align(str1_, str2_, cost_mat)
    alignment_ = result[2]
    alignment_.reverse()

    for (i, j) in alignment_:
        print(dec_map[i], dec_map[j])
    print(result[1])
