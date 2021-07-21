import leven
import pandas as pd
import numpy as np
from itertools import combinations

# load the transcription data
trs = pd.read_csv('transcriptions.tsv', sep='\t')

# read in the PMI distance matrix
enc_dict, dec_dict, dist_mat = leven.init_cost_matrix_weighted(
    'sounddistances_pmi.txt')
# note that the indices for characters are the keys in enc_dict
# as well as the values in dec_dict

# let's also rename the first column
trs = trs.rename(columns={'Unnamed: 0': 'location'})

# remember that transcriptions must be encoded as the indices of the matrix
trs_encoded = trs

for col in trs_encoded.iloc[:, 1:]:  # we skip the first column with locations
    trs_encoded[col] = [np.array([enc_dict[char] for char in trs]) if not pd.isna(
        trs) else trs for trs in trs_encoded[col]]

# now we can compute with the transcriptions
# let's compute the distance between all transcriptions of the first word

distances = []
alignments = []

for i, j in combinations(trs_encoded.iloc[:, 1], 2):
    # ignore missing values
    if not isinstance(i, float) and not isinstance(j, float):
        result = leven.leven_compute_align(i, j, dist_mat)

        # distance is the first value, the alignment in reverse is the third
        distances.append(result[0])
        alignments.append(result[2])

print(distances)
print(alignments)
