import leven
from itertools import combinations
import subprocess
import numpy as np
from rpy2.robjects.packages import importr
import rpy2.robjects as robjects
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
import os
from math import factorial
from time import time
from numba import njit


def bash(cmd):
    process = subprocess.run(cmd,
                             shell=True,
                             check=True,
                             universal_newlines=True,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT)
    return process


def get_align_from_sam(bash_output):
    alignments = set()

    split_ = bash_output.split('\n')
    for i, str_ in enumerate(split_):
        cur_align = []
        if len(str_) > 0:
            if str_[0] == '*':
                for segment in zip(str_[1:].strip(), split_[i + 1].strip()):
                    if not segment[0] == ' ':
                        cur_align.append(segment)
                alignments.add(tuple(cur_align))
    return alignments


def get_align_from_sam_3d(bash_output):
    alignments = set()

    split_ = bash_output.split('\n')
    for i, str_ in enumerate(split_):
        cur_align = []
        if len(str_) > 0:
            if str_[0] == '*':
                for segment in zip(str_[1:].strip(),
                                   split_[i + 1].strip(),
                                   split_[i + 2].strip()):
                    if not segment[0] == ' ':
                        cur_align.append(segment)
                alignments.add(tuple(cur_align))
    return alignments


def get_align_from_ipa(enc_align, decode_map):
    alignments = set()
    split_ = enc_align

    alignment = []
    for segment in split_:
        if segment == (-1, -1):
            alignment.reverse()
            alignments.add(tuple(alignment))
            alignment.clear()
        else:
            alignment.append(
                (decode_map[segment[0]], decode_map[segment[1]]))
    return alignments


def get_align_from_ipa_3d(enc_align, decode_map):
    alignments = set()
    split_ = enc_align

    alignment = []
    for segment in split_:
        if segment == (-1, -1, -1):
            alignment.reverse()
            alignments.add(tuple(alignment))
            alignment.clear()
        else:
            alignment.append(
                (decode_map[segment[0]],
                 decode_map[segment[1]],
                 decode_map[segment[2]]))
    return alignments


importr('ipa')
ipa = robjects.r['ipa']

with open('inventory.txt', 'r', encoding='utf8') as inv_file:
    char_inv = {line.strip().split('\t')[0]: line.strip().split('\t')[1]
                for line in inv_file.readlines()}
chars_ = list(char_inv.keys())
cost_mat = leven.init_cost_matrix(chars_, char_inv)
enc_map, dec_map = leven.generate_char_map(chars_)
dec_map[-1] = '-'

# init the numba compilation
leven.leven_compute_align(np.array([enc_map[char] for char in 'abc']),
                          np.array([enc_map[char] for char in 'abc']), cost_mat)

# init the numba compilation
leven.leven_3_dim(np.array([enc_map[char] for char in 'abc']),
                  np.array([enc_map[char] for char in 'abc']),
                  np.array([enc_map[char] for char in 'abc']), cost_mat)

# from encoded character-number to sampa symbol
dec_num_to_sam = {char: ipa(dec_map[char], to='xsampa')[0]
                  for char in dec_map}

# import data
df = pd.read_csv('trs_files/GTRP.txt', sep='\t')

for col in df.iloc[:, 1:]:
    df[col] = [leven.remove_accents(trs) for trs in df[col]]

unique_trs = list(
    {trs for col in df.iloc[:, 1:] for trs in df[col] if not pd.isna(trs)})
trs_trans = ipa(robjects.StrVector(unique_trs))
trans_ = {i: j for i, j in zip(unique_trs, trs_trans)}


def cycle(column, data_frame):
    trs_ = data_frame[column].dropna().unique()
    print(column, '\t', len(trs_))

    for i, j in combinations(trs_, 2):

        if 'ɫ' in i or 'ɫ' in j \
                or 'ɷ' in i or 'ɷ' in j \
                or 'ɼ' in i or 'ɼ' in j \
                or 'ɥ' in i or 'ɥ' in j \
                or 'ʍ' in i or 'ʍ' in j:
            continue
        else:
            # print(i, j)
            i_ = np.array([enc_map[char] for char in i])
            j_ = np.array([enc_map[char] for char in j])

            lev_result = leven.leven_compute_align(i_, j_, cost_mat)
            lev_output = get_align_from_ipa(lev_result[-1], dec_num_to_sam)

            i_sampa_ = trans_[i]
            j_sampa_ = trans_[j]
            if '\\' in i_sampa_ or '\\' in j_sampa_ or '&' in i_sampa_ \
                    or '&' in j_sampa_:
                continue
            else:
                bash_ = bash(
                    './3_dim_lev/levenshtein2 ' + i_sampa_ + ' ' + j_sampa_
                    + ' 3_dim_lev/data.cm0')
                wilbert_output = get_align_from_sam(bash_.stdout)

                diff = lev_output.difference(wilbert_output)

                if not len(diff) == 0:
                    print(i, j)
                    print(i_sampa_, j_sampa_)
                    print("ERROR")
                    exit()


def cycle_3d(column, data_frame):
    trs_ = data_frame[column].dropna().unique()
    print(column, '\t', len(trs_))

    for i, j, k in combinations(trs_, 3):
        if pd.isna(i) or pd.isna(j) or pd.isna(k)\
                or 'ɫ' in i or 'ɫ' in j or 'ɫ' in k \
                or 'ɷ' in i or 'ɷ' in j or 'ɷ' in k \
                or 'ɼ' in i or 'ɼ' in j or 'ɼ' in k \
                or 'ɥ' in i or 'ɥ' in j or 'ɥ' in k \
                or 'ʍ' in i or 'ʍ' in j or 'ʍ' in k:
            continue
        else:
            # print(column, ': ', i, j, k)
            i_ = np.array([enc_map[char] for char in i])
            j_ = np.array([enc_map[char] for char in j])
            k_ = np.array([enc_map[char] for char in k])

            lev_result = leven.leven_3_dim(i_, j_, k_, cost_mat)
            get_align_from_ipa_3d(lev_result[1], dec_num_to_sam)

            lev_output = get_align_from_ipa(lev_result[-1], dec_num_to_sam)
            i_sampa_ = trans_[i]
            j_sampa_ = trans_[j]
            k_sampa_ = trans_[k]
            if '\\' in i_sampa_ or '\\' in j_sampa_ or '\\' in k_sampa_ or \
                    '&' in i_sampa_ or '&' in j_sampa_ or '&' in k_sampa_:
                continue
            else:
                bash_ = bash(
                    './3_dim_lev/levenshtein3 '
                    + i_sampa_ + ' '
                    + j_sampa_ + ' '
                    + k_sampa_ + ' 3_dim_lev/data.cm0')
                wilbert_output = get_align_from_sam(bash_.stdout)

                diff = lev_output.difference(wilbert_output)

                if not len(diff) == 0:
                    print(i, j)
                    print(i_sampa_, j_sampa_)
                    print("ERROR")
                    print()
                    print(lev_output)
                    print(wilbert_output)
                    print()
                    exit()

# with ProcessPoolExecutor(
#         max_workers=len(os.sched_getaffinity(0))
#         # max_workers=1
# ) as executor:
#     results = {executor.submit(cycle, col, df) for col in df.columns[1:]}


with ProcessPoolExecutor(
        max_workers=len(os.sched_getaffinity(0)) - 2
        # max_workers=1
) as executor:
    try:
        results = {executor.submit(cycle_3d, col, df) for col in df.columns[1:]}
    except:
        print("ERROR")

# for col in df.columns[1:2]:
#     start = time()
#     cycle(col, df)
#     end = time()
#     print(col, '\t', end - start)
#     print()

# sum_fac = sum([cycle_3d(col, df) for col in df.columns[1:]])
# print(sum_fac)

# sum_fac = sum([cycle(col, df) for col in df.columns[1:]])
# print(sum_fac)

# s1 = 'ʔɒrde'
# s2 = 'ɛət'
# s3 = 'eɛrdə'

# s1_sam = trans_[s1]
# s2_sam = trans_[s2]
# s3_sam = trans_[s3]

# s1 = np.array([enc_map[char] for char in s1])
# s2 = np.array([enc_map[char] for char in s2])
# s3 = np.array([enc_map[char] for char in s3])

# lev_result = leven.leven_3_dim(s1, s2, s3, cost_mat)[1]
# bash_ = bash('./3_dim_lev/levenshtein3 ' + s1_sam + ' ' + s2_sam + ' '
#              + s3_sam + ' 3_dim_lev/data.cm0 1')

# print(get_align_from_ipa_3d(lev_result, dec_num_to_sam))
# print(get_align_from_sam_3d(bash_.stdout))
# print(bash_.stdout)
