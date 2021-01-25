import logging
from itertools import permutations, combinations
import numpy as np
import pandas as pd
from numba import njit
import sys

sys.path.insert(0, '..')
import leven

# logging.basicConfig(filename='computation.log', filemode='w', format='%(name)s - %(levelname)s - %(message)s',
#                     level=logging.DEBUG)
logging.basicConfig(filename='computation.log', filemode='w', format='%(name)s - %(levelname)s - %(message)s')


# np.set_printoptions(threshold=sys.maxsize)


def compute_all_cols(df, list_of_chars, trs_enc_map, cost_matrix):
    count_matrix = np.zeros((len(list_of_chars), len(list_of_chars)), dtype='int')
    # count = 0

    for col in df:
        cache = {}
        # count += 1
        for target_trs, compare_trs in combinations(df[col], 2):
            # only include strings, not NAs
            if target_trs in trs_enc_map and compare_trs in trs_enc_map:
                # use cached result if possible
                if (target_trs, compare_trs) not in cache:
                    result = leven.leven_compute_align(trs_enc_map[target_trs],
                                                       trs_enc_map[compare_trs], cost_matrix)
                    matrix_count(result[2], count_matrix)
                    cache[(target_trs, compare_trs)] = result[2]
                else:
                    matrix_count(cache[(target_trs, compare_trs)], count_matrix)
        # print(count)
    return count_matrix


@njit(cache=True)
def matrix_count(segment_list, count_mat):
    """ Update a segment-segment count matrix based on the indices stored in alignments"""
    for update_segment in segment_list:
        count_mat[update_segment] += 1
    return count_mat


@njit(cache=True)
def pmi(i, j, x_count, y_count, corpus_size, count_mat):
    # P(X, Y)
    # co_count = count_mat[i, j] + count_mat[j, i]
    co_count = count_mat[i, j]
    co_prob = co_count / corpus_size

    # P(X)
    x_prob = x_count / corpus_size

    # P(Y)
    y_prob = y_count / corpus_size

    if co_count == 0:
        return np.log2((0.1 ** 200) / (x_prob * y_prob))
    else:
        return np.log2(co_prob / (x_prob * y_prob))


def compute_pmi(count_matrix, list_of_chars):
    # matching operations are not taken into consideration for the corpus size
    np.fill_diagonal(count_matrix, 0)

    # compute row-wise, column-wise, and total sums
    total_sums, col_sums, row_sums = count_matrix.sum(), count_matrix.sum(axis=0), count_matrix.sum(axis=1)

    char_sums = {char: (col_sums[index] + row_sums[index]) for index, char in enumerate(list_of_chars)}
    pmi_mat = np.zeros(count_matrix.shape)

    for (i, j) in np.ndindex(count_matrix.shape):
        pmi_mat[i, j] = pmi(i, j,
                            char_sums[list_of_chars[i]],
                            char_sums[list_of_chars[j]],
                            total_sums, count_matrix)
    max_pmi_val = np.max(pmi_mat)

    # print(np.round(pmi_mat[-1]))

    for (i, j) in np.ndindex(pmi_mat.shape):
        pmi_mat[i, j] = 0 - pmi_mat[i, j] + max_pmi_val
    segment_vals = np.divide(pmi_mat, np.max(pmi_mat))

    # matches are set to a segment distance of 0
    np.fill_diagonal(segment_vals, 0)
    return segment_vals


def iteration_procedure(transcriptions, transcriptions_mapping, cost_matrix, iteration):
    # updated_results = get_dist(transcriptions, transcriptions_mapping, cost_matrix)
    updated_counts = compute_all_cols(trs, char_inv, trs_map, cost_mat)
    updated_pmi_matrix = compute_pmi(updated_counts, char_inv)
    print(updated_pmi_matrix)

    np.savetxt('cost_mat_iter_' + str(iteration) + '.txt', cost_mat - updated_pmi_matrix, delimiter='\t', fmt='%.4f')

    updated_difference = np.sum((cost_matrix - updated_pmi_matrix) ** 2)
    return updated_pmi_matrix, updated_difference


if __name__ == '__main__':
    dat = pd.read_csv('../trs_files/RND_no_dia.txt', sep='\t')
    # dat = pd.read_csv('GTRP.txt', sep='\t')
    trs = dat.iloc[:, 1:]

    trs = trs.iloc[:, 0:20]

    ''' Get rid of diacritics '''
    trs = pd.DataFrame.from_dict({col: [leven.remove_accents(transcription)
                                        for transcription in trs[col]]
                                  for col in trs})

    with open('../inventory.txt', 'r', encoding='utf8') as file:
        symbol_spec = {line.strip().split('\t')[0]: line.strip().split('\t')[1]
                       for line in file.readlines()}

    ''' generate an initial cost matrix with linguistic constraints'''
    char_inv = sorted(list(leven.get_char_inv(trs)))
    char_inv.append("-")  # add the indel character
    cost_mat = leven.init_cost_matrix(char_inv, symbol_spec)
    # np.savetxt('test_cost_mat', cost_mat, delimiter='\t', fmt='%.2f')

    # enc_trs, dec_trs, cost_mat = leven.init_cost_matrix_weighted('../pmi/pmi_wieling_large.tsv')

    # generate mappings to numpy arrays
    char_enc, char_dec = leven.generate_char_map(char_inv)
    trs_map = leven.encode_df(trs, char_enc)

    """ Levenshtein distance computations """
    # start = time.time()
    # results = get_dist(trs, trs_map, cost_mat)

    """ Count segment-segment occurrences """
    # counts = count_from_dist(results)
    # end = time.time()
    # print(end - start)
    """ Obtain new cost matrix (counts -> PMI -> segment distances) """
    # pmi_vals = compute_pmi(counts, char_inv)

    # np.savetxt('cost_mat_iter_2.txt', pmi_vals, delimiter='\t', fmt='%.8f')

    iteration, difference = 0, 999

    # # while difference != 0:
    while iteration < 10:
        iteration += 1
        cost_mat, difference = iteration_procedure(trs, trs_map, cost_mat, iteration)
        # np.savetxt('cost_mat_iter_' + str(iteration) + '.txt', cost_mat, delimiter='\t', fmt='%.4f')

        print('Iteration', iteration)
        # print(cost_mat)
        print(difference)
        print()

    # print(char_inv)

    # np.savetxt('converged_output.txt', cost_mat, fmt='%.8f', delimiter='\t')
    # cost_mat = np.loadtxt('converged_output.txt', delimiter='\t')

    ''' Testing alignments '''
    #
    #
    # char_inv = [char for char in symbol_spec]
    # cost_mat = leven.init_cost_matrix(char_inv, symbol_spec)
    # enc_trs, dec_trs = leven.generate_char_map(char_inv)

    # s1 = 'abdf'
    # s2 = 'abef'

    # s1 = np.array([enc_trs[char] for char in s1])
    # s2 = np.array([enc_trs[char] for char in s2])
    # print(s1)
    # print(s2)

    # dist = leven.leven_compute_align(s1, s2, cost_mat)
    # print(dist[1])

    # alignment = dist[2]
    # # print(alignment)
    # alignment.reverse()

    # # print(cost_mat[enc_trs['d'], -1])

    # for segment in alignment:
    #     # print(segment)
    #     if segment[0] == -1:
    #         print("- <->", dec_trs[segment[1]])
    #     elif segment[1] == -1:
    #         print(dec_trs[segment[0]], "<-> -")
    #     else:
    #         print(dec_trs[segment[0]], '<->', dec_trs[segment[1]])
