import leven
from itertools import combinations
import subprocess
import numpy as np
from rpy2.robjects.packages import importr
import rpy2.robjects as robjects
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
import os


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


importr('ipa')
ipa = robjects.r['ipa']

with open('inventory.txt', 'r', encoding='utf8') as inv_file:
    char_inv = {line.strip().split('\t')[0]: line.strip().split('\t')[1]
                for line in inv_file.readlines()}
chars_ = list(char_inv.keys())
cost_mat = leven.init_cost_matrix(chars_, char_inv)
enc_map, dec_map = leven.generate_char_map(chars_)
dec_map[-1] = '-'

dec_num_to_sam = {char: ipa(dec_map[char], to='xsampa')[0]
                  for char in dec_map}

# import data and clean it
df = pd.read_csv('trs_files/GTRP.txt', sep='\t')
for col in df.iloc[:, 1:]:
    df[col] = [leven.remove_accents(trs) for trs in df[col]]


def cycle(col, data_frame):
    unq = list(data_frame[col].dropna().unique())

    for i, j in combinations(unq, 2):
        # print(i, j)
        if 'ɫ' in i or 'ɫ' in j or 'ɷ' in i or 'ɷ' in j or 'ɼ' in i or \
                'ɼ' in j or 'ɥ' in i or 'ɥ' in j or 'ʍ' in i or 'ʍ' in j:
            pass
        else:
            str1 = np.array([enc_map[char] for char in i])
            str2 = np.array([enc_map[char] for char in j])
            lev_result = leven.leven_compute_align(str1, str2, cost_mat)
            lev_output = get_align_from_ipa(lev_result[-1], dec_num_to_sam)

            sam1 = ''.join([dec_num_to_sam[char] for char in str1])
            sam2 = ''.join([dec_num_to_sam[char] for char in str2])

            if '\\' in sam1 or '\\' in sam2 or '&' in sam1 or '&' in sam2:
                pass
            else:
                bash_ = bash('./3_dim_lev/levenshtein2 ' + sam1 +
                             ' ' + sam2 + ' 3_dim_lev/data.cm0')
                wilbert_output = get_align_from_sam(bash_.stdout)

                diff = lev_output.difference(wilbert_output)
                if not len(diff) == 0:
                    print(i, j)
                    print(sam1, sam2)
                    print("ERROR")
                    exit()
    print(col)


with ProcessPoolExecutor(
        max_workers=len(os.sched_getaffinity(0))
        # max_workers=1
) as executor:
    results = {executor.submit(cycle, col, df) for col in df.columns[1:4]}
