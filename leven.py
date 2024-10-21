""" compute VC sensitive Levenshtein distance cost matrix
    retrace cost matrix to find VC sensitive alignment"""
import functools
import itertools
from sys import maxsize
import numpy as np
import pandas as pd
from numba import int64, float64
from numba import njit
from numba.core.types.containers import UniTuple
from numba.typed import List
from unidecode import unidecode

list_type = UniTuple(int64, 2)
list_type_3d = UniTuple(int64, 3)


def get_char_inv(df, drop_col=None):
    """get the character inventory of a DataFrame;
    possible to remove a location column with the drop_col argument"""
    if drop_col:
        return {char for trs in df.drop(columns=drop_col).stack() for char in trs}
    else:
        return {char for trs in df.stack() for char in trs}


def remove_accents(input_str: str):
    """removes known diacritics from unicode type strings;
    note that currently the '/' and spaces are also removed,
    although they are used to denote multiple transcriptions
    in Gabmap style data"""
    if isinstance(input_str, str):
        # if two transcriptions are given (in Gabmap format), take the first
        if "/" in input_str:
            input_str = input_str.split("/")[0]

        remove_chars_manual = {
            " ",
            "_",
            "/",
            ":",
            "´",
            "ʰ",
            "ʲ",
            "ʷ",
            "ˀ",
            "ˈ",
            "ː",
            "ˑ",
            "˞",
            "ˠ",
            "ˡ",
            "ˢ",
            "ˤ",
            "˳",
            "˺",
            "̂",
            "̃",
            "̆",
            "̌",
            "̍",
            "̑",
            "̘",
            "̙",
            "̚",
            "̜",
            "̝",
            "̞",
            "̟",
            "̠",
            "̤",
            "̥",
            "̩",
            "̪",
            "̬",
            "̮",
            "̯",
            "̰",
            "̱",
            "̹",
            "̺",
            "̻",
            "͂",
            "͆",
            "͉",
            "͎",
            "◌̭",
            "͡",
            "ᵊ",
            "̴",
            "ʲ",
        }

        # manually remove chars
        final_form = "".join(
            char for char in input_str if char not in remove_chars_manual
        )
        return final_form
    else:
        return np.nan


def string2ascii(input_str: str):
    """removes diacritics from unicode strings, leaving only ASCII symbols;
    note that currently the '/' and spaces are also removed,
    although they are used to denote multiple transcriptions
    in Gabmap style data"""
    if isinstance(input_str, str):
        # if two transcriptions are given (in Gabmap format), take the first
        if "/" in input_str:
            input_str = input_str.split("/")[0]
        return unidecode(input_str)
    else:
        return np.nan


def generate_trs_map(list_of_trs, char_map):
    """generates NumPy number encodings given a list"""
    return {
        trs: np.array([char_map[char] for char in trs])
        for trs in list_of_trs
        if isinstance(trs, str)
    }


def encode_df(df, char_map):
    """generates NumPy number encodings given a dataframe"""
    return {
        trs: np.array([char_map[char] for char in trs])
        for trs in df.stack()
        if isinstance(trs, str)
    }


def generate_char_map(list_of_chars):
    """generate a character mapping from list of characters"""
    enc_dict = {char: index for index, char in enumerate(list_of_chars)}
    dec_dict = dict(zip(enc_dict.values(), enc_dict.keys()))
    return enc_dict, dec_dict


def encode_column_trs_array(list_of_trs, char_map):
    """encode a Pandas DataFrame column of transcriptions to NumPy arrays"""
    encoded_list = []
    for trs in list_of_trs:
        if isinstance(trs, float):
            encoded_list.append(np.array([-1]))
        else:
            encoded_list.append(np.array([char_map[char] for char in trs]))
    return encoded_list


def encode_column_trs_array_native(list_of_trs, char_map):
    """encode a Pandas DataFrame column of transcriptions to NumPy arrays"""
    encoded_list = []
    for trs in list_of_trs:
        if isinstance(trs, float):
            encoded_list.append(np.nan)
        else:
            encoded_list.append(np.array([char_map[char] for char in trs]))
    return encoded_list


def encode_column_trs(list_of_trs, char_map):
    """encode a Pandas DataFrame column of transcriptions to NumPy arrays"""
    encoded_list = []
    for trs in list_of_trs:
        if isinstance(trs, float):
            encoded_list.append(np.nan)
        else:
            encoded_list.append(np.array([char_map[char] for char in trs]))
    return encoded_list


def init_cost_matrix_semivowels(list_of_chars, symbol_specification_dict):
    """generates a NumPy segment-segment matrix  **without** an input file
    semivowels are allowed to match both vowels and consonants"""
    inv_size = len(list_of_chars)
    cost_matrix = np.zeros((inv_size, inv_size))

    for i, j in np.ndindex(inv_size, inv_size):
        # character inventory indices are used for filling in the cost matrix
        c1 = list_of_chars[i]
        c2 = list_of_chars[j]

        if c1 == c2:
            # match
            cost_matrix[i, j] = 0
        elif c1 == "-" or c2 == "-":
            # indel
            if c1 == "-":
                cost_matrix[-1, j] = 1
            if c2 == "-":
                cost_matrix[i, -1] = 1
        # elif c1 in {'j', 'w', 'i', 'u'} or c2 in {'j', 'w', 'i', 'u'}:
        #     cost_matrix[i, j] = 1
        elif c1 == "ə" or c2 == "ə":
            if c1 in {"m", "l", "n", "r", "ŋ", "j", "w"} or c2 in {
                "m",
                "l",
                "n",
                "r",
                "ŋ",
                "j",
                "w",
            }:
                cost_matrix[i, j] = 1
            else:
                if (
                    c1 != c2
                    and symbol_specification_dict[c1] == symbol_specification_dict[c2]
                ):
                    # substitution
                    cost_matrix[i, j] = 1
                elif (
                    c1 != c2
                    and symbol_specification_dict[c1] != symbol_specification_dict[c2]
                ):
                    # V-C substitution
                    cost_matrix[i, j] = 100
        elif (
            c1 != c2 and symbol_specification_dict[c1] == symbol_specification_dict[c2]
        ):
            # substitution
            cost_matrix[i, j] = 1
        elif (
            c1 != c2 and symbol_specification_dict[c1] != symbol_specification_dict[c2]
        ):
            # V-C substitution
            cost_matrix[i, j] = 100
    return cost_matrix


def init_cost_matrix(list_of_chars, symbol_specification_dict):
    """generates a NumPy segment-segment matrix  **without** an input file"""
    inv_size = len(list_of_chars)
    cost_matrix = np.zeros((inv_size, inv_size))

    for i, j in np.ndindex(inv_size, inv_size):
        # character inventory indices are used for filling in the cost matrix
        c1 = list_of_chars[i]
        c2 = list_of_chars[j]

        if c1 == c2:
            # match
            cost_matrix[i, j] = 0
        elif c1 == "-" or c2 == "-":
            # indel
            if c1 == "-":
                cost_matrix[-1, j] = 1
            if c2 == "-":
                cost_matrix[i, -1] = 1
        elif (
            c1 != c2 and symbol_specification_dict[c1] == symbol_specification_dict[c2]
        ):
            # substitution
            cost_matrix[i, j] = 1
        elif (
            c1 != c2 and symbol_specification_dict[c1] != symbol_specification_dict[c2]
        ):
            # V-C substitution
            cost_matrix[i, j] = 100
    return cost_matrix


def init_cost_matrix_weighted(pmi_file):
    """generates a NumPy segment-segment matrix  **with** an input file"""
    with open(pmi_file, "r") as file:
        char_list = file.readlines()[0].strip().split("\t")
        enc_dict, dec_dict = generate_char_map(char_list)

    dist_mat = pd.read_csv(pmi_file, sep="\t").drop(columns="Unnamed: 0").to_numpy()
    return enc_dict, dec_dict, dist_mat


@njit(cache=True)
def weight_(c1_idx, c2_idx, c3_idx, cost_matrix):
    """compute the sum of all pairwise weights"""
    return (
        cost_matrix[c1_idx, c2_idx]
        + cost_matrix[c1_idx, c3_idx]
        + cost_matrix[c2_idx, c3_idx]
    )


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
    """check if current option is the final one"""
    # return (node[0], node[1]) == (valid_vals[-2], valid_vals[-1])
    # return node == valid_vals[-2::]

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
def check_op(op, op_vec, w1_array, w2_array, trace_list, align_, trace):
    if op == 0:
        align_.append((w1_array[op_vec[2]], -1))
        if trace:
            trace_list.insert(0, (op_vec[2], op_vec[3]))
    if op == 1:
        align_.append((-1, w2_array[op_vec[5]]))
        if trace:
            trace_list.insert(0, (op_vec[4], op_vec[5]))
    if op == 2:
        align_.append((w1_array[op_vec[6]], w2_array[op_vec[7]]))
        if trace:
            trace_list.insert(0, (op_vec[6], op_vec[7]))


@njit(cache=True)
def backtrack(a, array, w1_array, w2_array, trace_list, align_, alignment):
    """backtrack the path from the pointer array without strings;
    faster for PMI matrix computation,
    although strings still need to be decoded for human legibility"""
    # construct reshaped array of only source node coordinates
    iter_ = a.reshape(4, 2)[1:]

    # compute the number of options (note: value is doubled)
    options = np.sum(np.not_equal(iter_, -1))

    # view only valid coordinates for index checking (see below)
    valid_ = a[a != -1]

    for i, k in enumerate(iter_):
        # value is [-1, -1], so continue iteration
        if np.max(k) == -1:
            continue
        # value is [0, 0], so we are at the beginning
        if np.max(k) == 0:
            check_op(i, a, w1_array, w2_array, trace_list, align_, False)
            align_.append((-1, -1))
            add_alignment(align_, alignment)

            if len(trace_list) == 0:
                return 0
            else:
                align_.clear()
                return backtrack(
                    array[-1], array, w1_array, w2_array, trace_list, align_, alignment
                )
        else:
            check_op(i, a, w1_array, w2_array, trace_list, align_, False)
            return backtrack(
                array[(array[:, 0] == k[0]) & (array[:, 1] == k[1])].flatten(),
                array,
                w1_array,
                w2_array,
                trace_list,
                align_,
                alignment,
            )


@njit(cache=True)
def decompose(alignment, cost_mat):
    var1_var2 = 0

    for segment in alignment:
        if segment == (-1, -1):
            break  # only first alignment, as the distance is always the same
        else:
            var1_var2 += cost_mat[segment[0], segment[1]]

    max_length = 0
    length = 0
    for segment in alignment:
        if segment == (-1, -1):
            if length > max_length:  # take the longest alignment length
                max_length = length
            length = 0
        else:
            length += 1

    return var1_var2, max_length


@njit(cache=True)
def leven_compute_align(w1_idx, w2_idx, cost_matrix):
    """compute the dynamic programming tableau;
    track 'pointers' from cell to cell"""
    row_nr = 0

    # init tableau
    # tabl = np.zeros((w1_idx.size + 1, w2_idx.size + 1))
    tabl = np.empty((w1_idx.size + 1, w2_idx.size + 1))
    tabl[0, 0] = 0

    # initiate the pointer matrix with lengths corresponding to input strings
    path_directions = np.full(((tabl.size - 1), 8), -1)

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
    for i, j in np.ndindex(w1_idx.size + 1, w2_idx.size + 1):
        # ignore the first column and row
        if i == 0 or j == 0:
            continue
        path_directions[row_nr][0] = i
        path_directions[row_nr][1] = j

        # compute prior costs + added operation costs
        # going down: ins s1
        ins_s1 = tabl[i - 1, j] + cost_matrix[-1, w1_idx[i - 1]]

        # going right: ins s2
        ins_s2 = tabl[i, j - 1] + cost_matrix[-1, w2_idx[j - 1]]

        # going down-right
        sub = tabl[i - 1, j - 1] + cost_matrix[w1_idx[i - 1], w2_idx[j - 1]]
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

        tabl[i, j] = min_val
        row_nr += 1

    """ 3 lists are required: the final set of alignments, one temporary list
        for each alignment, and one list for tracing visited crossings"""
    alignment = List.empty_list(item_type=list_type)
    current_ = List.empty_list(item_type=list_type)
    trace = List.empty_list(item_type=list_type)

    backtrack(
        path_directions[-1], path_directions, w1_idx, w2_idx, trace, current_, alignment
    )
    dists = decompose(alignment, cost_matrix)
    return tabl[-1, -1], tabl, alignment, dists, path_directions


@njit(cache=True)
def min_(o1, o2, o3, o4, o5, o6, o7):
    """returns the minimum value from the given operation costs"""
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
def check_first_3d(node, list_):
    if len(list_) == 0:
        return False
    else:
        if (node[0], node[1], node[2]) == list_[0]:
            return True
        else:
            return False


@njit(cache=True)
def check_depth_3d(node, trace_list):
    if len(trace_list) > 0:
        if (
            trace_list[0][0] < node[0]
            or trace_list[0][1] < node[1]
            or trace_list[0][2] < node[2]
        ):
            return True
        else:
            return False
    else:
        return False


@njit(cache=True)
def check_final_3d(node, valid_vals):
    """check if current option is the final one"""
    if (node[0], node[1], node[2]) == (valid_vals[-3], valid_vals[-2], valid_vals[-1]):
        return True
    else:
        return False


@njit(cache=True)
def member_3d_(node, trace_list):
    is_visited = False
    for check_node in trace_list:
        if (node[0], node[1], node[2]) == check_node:
            is_visited = True
    return is_visited


@njit(cache=True)
def check_op_3d(op, op_vec, w1_array, w2_array, w3_array, trace_list, align_, trace):
    if op == 0:
        align_.append((w1_array[op_vec[5]], -1, -1))
        if trace:
            trace_list.insert(0, (op_vec[3], op_vec[4], op_vec[5]))
    if op == 1:
        align_.append((-1, w2_array[op_vec[7]], -1))
        if trace:
            trace_list.insert(0, (op_vec[6], op_vec[7], op_vec[8]))
    if op == 2:
        align_.append((-1, -1, w3_array[op_vec[9]]))
        if trace:
            trace_list.insert(0, (op_vec[9], op_vec[10], op_vec[11]))
    if op == 3:
        align_.append((-1, w2_array[op_vec[13]], w3_array[op_vec[12]]))
        if trace:
            trace_list.insert(0, (op_vec[12], op_vec[13], op_vec[14]))
    if op == 4:
        align_.append((w1_array[op_vec[17]], -1, w3_array[op_vec[15]]))
        if trace:
            trace_list.insert(0, (op_vec[15], op_vec[16], op_vec[17]))
    if op == 5:
        align_.append((w1_array[op_vec[20]], w2_array[op_vec[19]], -1))
        if trace:
            trace_list.insert(0, (op_vec[18], op_vec[19], op_vec[20]))
    if op == 6:
        align_.append(
            (w1_array[op_vec[23]], w2_array[op_vec[22]], w3_array[op_vec[21]])
        )
        if trace:
            trace_list.insert(0, (op_vec[21], op_vec[22], op_vec[23]))


@njit(cache=True)
def decompose_3d(alignment, cost_mat, std_comp=False):
    if std_comp:
        var1_var2 = maxsize
        var2_var3 = maxsize
        var1_var3 = maxsize

        var1_var2_ = 0
        var2_var3_ = 0
        var1_var3_ = 0

        for triplet in alignment:
            if triplet == (-1, -1, -1):
                # if current 2-3 or 1-3 dist is lower, replace vals
                if var2_var3_ < var2_var3 or var1_var3_ < var1_var3:
                    var2_var3 = var2_var3_
                    var1_var3 = var1_var3_
                    var1_var2 = var1_var2_
                var1_var2_ = 0
                var2_var3_ = 0
                var1_var3_ = 0
            else:
                var1_var2_ += cost_mat[triplet[0], triplet[1]]
                var2_var3_ += cost_mat[triplet[1], triplet[2]]
                var1_var3_ += cost_mat[triplet[0], triplet[2]]
    else:
        var1_var2 = 0
        var2_var3 = 0
        var1_var3 = 0

        for triplet in alignment:
            if triplet == (-1, -1, -1):
                break  # only first alignment, as the distance is always the
                # same
            else:
                var1_var2 += cost_mat[triplet[0], triplet[1]]
                var2_var3 += cost_mat[triplet[1], triplet[2]]
                var1_var3 += cost_mat[triplet[0], triplet[2]]

    max_length = 0
    length = 0
    for triplet in alignment:
        if triplet == (-1, -1, -1):
            if length > max_length:  # take the longest alignment length
                max_length = length
            length = 0
        else:
            length += 1

    return var1_var2, var2_var3, var1_var3, max_length


@njit(cache=True)
def decompose_3d_count(alignment, cost_mat):
    var1_var2 = maxsize
    var2_var3 = maxsize
    var1_var3 = maxsize

    var1_var2_ = 0
    var2_var3_ = 0
    var1_var3_ = 0

    count_c = 0
    count_d = 0
    count_n = 0

    count_c_ = 0
    count_d_ = 0
    count_n_ = 0
    for triplet in alignment:
        if triplet == (-1, -1, -1):
            # if current 2-3 or 1-3 dist is lower, replace vals
            if var2_var3_ < var2_var3 or var1_var3_ < var1_var3:
                var2_var3 = var2_var3_
                var1_var3 = var1_var3_
                var1_var2 = var1_var2_

                count_c = count_c_
                count_d = count_d_
                count_n = count_n_
            var1_var2_ = 0
            var2_var3_ = 0
            var1_var3_ = 0

            count_c_ = 0
            count_d_ = 0
            count_n_ = 0
        else:
            old_new = cost_mat[triplet[0], triplet[1]]
            new_std = cost_mat[triplet[1], triplet[2]]
            old_std = cost_mat[triplet[0], triplet[2]]

            convergence = old_std - new_std

            var1_var2_ += old_new
            var2_var3_ += new_std
            var1_var3_ += old_std

            if convergence == 0:
                count_n_ += 1
                # count_n_ += abs(old_new)
            elif convergence > 0:
                # count_c_ += 1
                count_c_ += abs(convergence)
            elif convergence < 0:
                # count_d_ += 1
                count_d_ += abs(convergence)
            else:
                break

    max_length = 0
    length = 0
    for triplet in alignment:
        if triplet == (-1, -1, -1):
            if length > max_length:  # take the longest alignment length
                max_length = length
            length = 0
        else:
            length += 1
    return var1_var2, var2_var3, var1_var3, max_length, count_n, count_d, count_c


@njit(cache=True)
def decompose_3d_count_neutral_char(alignment, cost_mat, neutral_char=None):
    var1_var2 = maxsize
    var2_var3 = maxsize
    var1_var3 = maxsize

    var1_var2_ = 0
    var2_var3_ = 0
    var1_var3_ = 0

    count_c = 0
    count_d = 0
    count_n = 0

    count_c_ = 0
    count_d_ = 0
    count_n_ = 0
    for triplet in alignment:
        if triplet == (-1, -1, -1):
            # if current 2-3 or 1-3 dist is lower, replace vals
            if var2_var3_ < var2_var3 or var1_var3_ < var1_var3:
                var2_var3 = var2_var3_
                var1_var3 = var1_var3_
                var1_var2 = var1_var2_

                count_c = count_c_
                count_d = count_d_
                count_n = count_n_
            var1_var2_ = 0
            var2_var3_ = 0
            var1_var3_ = 0

            count_c_ = 0
            count_d_ = 0
            count_n_ = 0
        else:
            if (
                (triplet[0] == neutral_char)
                or (triplet[1] == neutral_char)
                or (triplet[2] == neutral_char)
            ):
                old_new = 0
                new_std = 0
                old_std = 0

                convergence = old_std - new_std

                var1_var2_ += old_new
                var2_var3_ += new_std
                var1_var3_ += old_std

                if convergence == 0:
                    count_n_ += 1
                    # count_n_ += abs(old_new)
                elif convergence > 0:
                    # count_c_ += 1
                    count_c_ += abs(convergence)
                elif convergence < 0:
                    # count_d_ += 1
                    count_d_ += abs(convergence)
                else:
                    break
            else:
                old_new = cost_mat[triplet[0], triplet[1]]
                new_std = cost_mat[triplet[1], triplet[2]]
                old_std = cost_mat[triplet[0], triplet[2]]

                convergence = old_std - new_std

                var1_var2_ += old_new
                var2_var3_ += new_std
                var1_var3_ += old_std

                if convergence == 0:
                    count_n_ += 1
                    # count_n_ += abs(old_new)
                elif convergence > 0:
                    # count_c_ += 1
                    count_c_ += abs(convergence)
                elif convergence < 0:
                    # count_d_ += 1
                    count_d_ += abs(convergence)
                else:
                    break

    max_length = 0
    length = 0
    for triplet in alignment:
        if triplet == (-1, -1, -1):
            if length > max_length:  # take the longest alignment length
                max_length = length
            length = 0
        else:
            length += 1
    return var1_var2, var2_var3, var1_var3, max_length, count_n, count_d, count_c


@njit(cache=True)
def leven_3_dim(w1_idx, w2_idx, w3_idx, cost_matrix, std_comp=False):
    """compute the dynamic programming tableau;
    track 'pointers' from cell to cell"""
    row_nr = 0

    # tableau initialization
    tabl = np.zeros((w3_idx.size + 1, w2_idx.size + 1, w1_idx.size + 1))

    # initiate the pointer matrix with lengths corresponding to input strings
    path_directions = np.full(((tabl.size - 1), 24), -1)

    # begin iteration
    for (i, j, k), val in np.ndenumerate(tabl):
        path_directions[row_nr][0] = i
        path_directions[row_nr][1] = j
        path_directions[row_nr][2] = k

        # init the operation costs to a large number
        i1 = maxsize
        i2 = maxsize
        i3 = maxsize
        i12 = maxsize
        i13 = maxsize
        i23 = maxsize
        i123 = maxsize

        # skip the origin
        if max(i, j, k) == 0:
            continue

        # determine costs for each operation
        if k > 0:
            # ins s1
            i1 = tabl[i, j, k - 1] + weight_(w1_idx[k - 1], -1, -1, cost_matrix)
        if j > 0:
            # ins s2
            i2 = tabl[i, j - 1, k] + weight_(-1, w2_idx[j - 1], -1, cost_matrix)
        if i > 0:
            # ins s3
            i3 = tabl[i - 1, j, k] + weight_(-1, -1, w3_idx[i - 1], cost_matrix)
        if j > 0 and i > 0:
            # ins s2 + s3
            i23 = tabl[i - 1, j - 1, k] + weight_(
                -1, w2_idx[j - 1], w3_idx[i - 1], cost_matrix
            )

        if k > 0 and i > 0:
            # ins s1 + s3
            i13 = tabl[i - 1, j, k - 1] + weight_(
                w1_idx[k - 1], -1, w3_idx[i - 1], cost_matrix
            )
        if k > 0 and j > 0:
            # ins s1 + s2
            i12 = tabl[i, j - 1, k - 1] + weight_(
                w1_idx[k - 1], w2_idx[j - 1], -1, cost_matrix
            )
        if i > 0 and j > 0 and k > 0:
            # ins s1 + s2 + s3
            i123 = tabl[i - 1, j - 1, k - 1] + weight_(
                w1_idx[k - 1], w2_idx[j - 1], w3_idx[i - 1], cost_matrix
            )
        min_val = min_(i1, i2, i3, i23, i13, i12, i123)
        tabl[i, j, k] = min_val

        if i1 == min_val:
            path_directions[row_nr][3] = i
            path_directions[row_nr][4] = j
            path_directions[row_nr][5] = k - 1
        if i2 == min_val:
            path_directions[row_nr][6] = i
            path_directions[row_nr][7] = j - 1
            path_directions[row_nr][8] = k
        if i3 == min_val:
            path_directions[row_nr][9] = i - 1
            path_directions[row_nr][10] = j
            path_directions[row_nr][11] = k
        if i23 == min_val:
            path_directions[row_nr][12] = i - 1
            path_directions[row_nr][13] = j - 1
            path_directions[row_nr][14] = k
        if i13 == min_val:
            path_directions[row_nr][15] = i - 1
            path_directions[row_nr][16] = j
            path_directions[row_nr][17] = k - 1
        if i12 == min_val:
            path_directions[row_nr][18] = i
            path_directions[row_nr][19] = j - 1
            path_directions[row_nr][20] = k - 1
        if i123 == min_val:
            path_directions[row_nr][21] = i - 1
            path_directions[row_nr][22] = j - 1
            path_directions[row_nr][23] = k - 1
        row_nr += 1

    """ 3 lists are required: the final set of alignments, one temporary list
            for each alignment, and one list for tracing visited crossings"""
    alignment = List.empty_list(item_type=list_type_3d)
    current_ = List.empty_list(item_type=list_type_3d)
    trace = List.empty_list(item_type=list_type_3d)
    backtrack_3d(
        path_directions[-1],
        path_directions,
        w1_idx,
        w2_idx,
        w3_idx,
        trace,
        current_,
        alignment,
    )
    if std_comp:
        dists = decompose_3d(alignment, cost_matrix, std_comp=True)
    else:
        dists = decompose_3d(alignment, cost_matrix)
    return tabl, alignment, dists


@njit(cache=True)
def backtrack_3d(a, array, w1_array, w2_array, w3_array, trace_list, align_, alignment):
    """backtrack the path from the pointer array without strings;
    faster for PMI matrix computation,
    although strings still need to be decoded for human legibility"""
    # construct reshaped array of only source node coordinates
    iter_ = a.reshape(8, 3)[1:]

    # compute the number of options (note: value is tripled)
    options = np.sum(np.not_equal(iter_, -1))

    # view only valid coordinates for index checking (see below)
    valid_ = a[a != -1]

    for i, k in enumerate(iter_):
        # value is [-1, -1], so continue iteration
        if np.max(k) == -1:
            continue
        # value is [0, 0], so we are at the beginning
        if np.max(k) == 0:
            check_op_3d(i, a, w1_array, w2_array, w3_array, trace_list, align_, False)
            align_.append((-1, -1, -1))
            add_alignment(align_, alignment)
            align_.clear()
            if len(trace_list) == 0:
                return 0
            else:
                return backtrack_3d(
                    array[-1],
                    array,
                    w1_array,
                    w2_array,
                    w3_array,
                    trace_list,
                    align_,
                    alignment,
                )
        else:
            if options > 3:
                final_ = check_final_3d(k, valid_)

                if member_3d_(k, trace_list):
                    if check_first_3d(k, trace_list):
                        trace_list.pop(0)
                    else:
                        check_op_3d(
                            i,
                            a,
                            w1_array,
                            w2_array,
                            w3_array,
                            trace_list,
                            align_,
                            False,
                        )
                        return backtrack_3d(
                            array[
                                (array[:, 0] == k[0])
                                & (array[:, 1] == k[1])
                                & (array[:, 2] == k[2])
                            ].flatten(),
                            array,
                            w1_array,
                            w2_array,
                            w3_array,
                            trace_list,
                            align_,
                            alignment,
                        )
                else:
                    if check_depth_3d(k, trace_list):
                        if final_:
                            check_op_3d(
                                i,
                                a,
                                w1_array,
                                w2_array,
                                w3_array,
                                trace_list,
                                align_,
                                False,
                            )
                            return backtrack_3d(
                                array[
                                    (array[:, 0] == k[0])
                                    & (array[:, 1] == k[1])
                                    & (array[:, 2] == k[2])
                                ].flatten(),
                                array,
                                w1_array,
                                w2_array,
                                w3_array,
                                trace_list,
                                align_,
                                alignment,
                            )
                    else:
                        if final_:
                            check_op_3d(
                                i,
                                a,
                                w1_array,
                                w2_array,
                                w3_array,
                                trace_list,
                                align_,
                                False,
                            )
                            return backtrack_3d(
                                array[
                                    (array[:, 0] == k[0])
                                    & (array[:, 1] == k[1])
                                    & (array[:, 2] == k[2])
                                ].flatten(),
                                array,
                                w1_array,
                                w2_array,
                                w3_array,
                                trace_list,
                                align_,
                                alignment,
                            )
                        else:
                            check_op_3d(
                                i,
                                a,
                                w1_array,
                                w2_array,
                                w3_array,
                                trace_list,
                                align_,
                                True,
                            )
                            return backtrack_3d(
                                array[
                                    (array[:, 0] == k[0])
                                    & (array[:, 1] == k[1])
                                    & (array[:, 2] == k[2])
                                ].flatten(),
                                array,
                                w1_array,
                                w2_array,
                                w3_array,
                                trace_list,
                                align_,
                                alignment,
                            )
            else:
                check_op_3d(
                    i, a, w1_array, w2_array, w3_array, trace_list, align_, False
                )
                return backtrack_3d(
                    array[
                        (array[:, 0] == k[0])
                        & (array[:, 1] == k[1])
                        & (array[:, 2] == k[2])
                    ].flatten(),
                    array,
                    w1_array,
                    w2_array,
                    w3_array,
                    trace_list,
                    align_,
                    alignment,
                )


@njit(cache=True)
def leven_dist(w1_idx, w2_idx, cost_matrix):
    """compute the dynamic programming tableau"""
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
        for i, j in np.ndindex(w1_idx.size + 1, w2_idx.size + 1):
            # ignore the first column and row
            if i != 0 and j != 0:
                # going down: ins s1 / del s2
                ins_s1 = tabl[i - 1, j] + cost_matrix[-1, w1_idx[i - 1]]

                # going right: del s1 / ins s2
                ins_s2 = tabl[i, j - 1] + cost_matrix[-1, w2_idx[j - 1]]

                # going down-right
                sub = tabl[i - 1, j - 1] + cost_matrix[w1_idx[i - 1], w2_idx[j - 1]]

                # compute prior costs + added operation costs
                costs = np.array([ins_s1, ins_s2, sub])

                # update tableau
                tabl[i, j] = np.min(costs)
        return tabl[-1, -1], tabl


def normalize_pmi(pmi_matrix, threshold=0.7):
    # obtain largest allowed segment-segment distance
    max_val = np.max(pmi_matrix[pmi_matrix < threshold])
    min_val = np.min(pmi_matrix)

    # normalize the matrix
    pmi_matrix = (pmi_matrix - min_val) / (max_val - min_val)
    return pmi_matrix


def test_2d_ipa():
    """2 dimensional test: IPA-based"""
    # generate cost matrix
    cost_mat = init_cost_matrix(chars_, char_inv)

    # encoder/decoder maps between numeric and character representations
    enc_map, dec_map = generate_char_map(chars_)
    dec_map[-1] = "-"

    str1 = "abc"
    str2 = "abe"

    # encode strings into NumPy arrays
    str1_ = np.array([enc_map[char] for char in str1])
    str2_ = np.array([enc_map[char] for char in str2])

    # compute Levenshtein distance
    result = leven_compute_align(str1_, str2_, cost_mat)

    # first value is distance
    print(result[0])

    # alignment is third result parameter, in reverse order, delimited by - -> -
    alignment = result[-3]
    alignment.reverse()
    for i in alignment:
        print(
            " -> ".join([dec_map[idx] for idx in i if not (i[0] == -1 and i[1] == -1)])
        )


def test_3d_ipa():
    """3 dimensional test: IPA-based"""
    # generate cost matrix
    cost_mat = init_cost_matrix(chars_, char_inv)

    # encoder/decoder maps between numeric and character representations
    enc_map, dec_map = generate_char_map(chars_)
    dec_map[-1] = "-"

    str1 = "abc"
    str2 = "abe"
    str3 = "aaa"

    # encode strings into NumPy arrays
    str1_ = np.array([enc_map[char] for char in str1])
    str2_ = np.array([enc_map[char] for char in str2])
    str3_ = np.array([enc_map[char] for char in str3])

    # compute Levenshtein distance
    result = leven_3_dim(str1_, str2_, str3_, cost_mat)

    # last value is tuple of distances (1-2, 2-3, 1-3) and alignment length
    print(result[-1])

    alignment = result[-2]
    alignment.reverse()

    for i in alignment:
        print(
            " <-> ".join(
                [
                    dec_map[idx]
                    for idx in i
                    if not (i[0] == -1 and i[1] == -1 and i[2] == -1)
                ]
            )
        )


def encode_graphemes(input_str: str, symbol_list: list, encoder_map: dict):
    """in case of multi-character input inventory;
    graphemes should be specified in inventory.txt"""

    encoded_str = []
    current_grapheme = input_str[0]
    idx = 0

    while idx < len(input_str):
        # exit at the last character
        if idx == len(input_str) - 1:
            encoded_str.append(encoder_map[current_grapheme])
            idx += 1
        else:
            # check if the next character occurs with current
            if current_grapheme + input_str[idx + 1] in symbol_list:
                current_grapheme += input_str[idx + 1]
                idx += 1
            # append current character(s) and continue
            else:
                encoded_str.append(encoder_map[current_grapheme])
                idx += 1
                current_grapheme = input_str[idx]

    return np.array(encoded_str)


def test_2d_grapheme():
    """2 dimensional test: grapheme-based"""
    # generate cost matrix
    cost_mat = init_cost_matrix(chars_, char_inv)

    # encoder/decoder maps between numeric and character representations
    enc_map, dec_map = generate_char_map(chars_)
    dec_map[-1] = "-"

    str1 = "gemeinskop"
    str2 = "gemeenscheeppen"

    # encode strings into NumPy arrays
    str1_ = encode_graphemes(str1, chars_, enc_map)
    str2_ = encode_graphemes(str2, chars_, enc_map)

    # compute Levenshtein distance
    result = leven_compute_align(str1_, str2_, cost_mat)

    # first value is distance
    print(result[0])

    # alignment is third result parameter, in reverse order, delimited by -
    # -> -
    alignment = result[-3]
    alignment.reverse()
    for i in alignment:
        print(
            " -> ".join([dec_map[idx] for idx in i if not (i[0] == -1 and i[1] == -1)])
        )


if __name__ == "__main__":
    # read in character/grapheme inventory
    with open("inventory.txt", "r", encoding="utf8") as inv_file:
        char_inv = {
            line.strip().split("\t")[0]: line.strip().split("\t")[1]
            for line in inv_file.readlines()
        }
    chars_ = list(char_inv.keys())

    test_2d_ipa()
    # test_3d_ipa()
    # test_2d_grapheme()
    #
    # with open(
    #         '/home/raoul/DOWNLOADS'
    #         '/nl_sassisk_nedderlandske_lemmaparen_filterd.txt',
    #         'r', encoding='utf8') as textfile:
    #     data = [[line.strip().split('\t')[0], line.strip().split('\t')[1]]
    #             for line in textfile.readlines()]
    # print([i + " -> " + string2ascii(i) for i in list(itertools.chain(*data))])
