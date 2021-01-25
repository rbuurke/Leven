from pathlib import Path

import pandas as pd

df_list = []

for file in Path('./matrices/gtrp_diareg').glob('*.txt'):
    print(file)
    with open(file, 'r') as open_file:
        output = open_file.readline().split('\t')
        location_sample_size = int(str(file).split('_')[3])
        word_sample_size = int(str(file).split('_')[-1].split('.')[0])
        for num in output:
            df_list.append([location_sample_size, word_sample_size, num])

df = pd.DataFrame(df_list)
df.columns = ['loc_ssize', 'word_ssize', 'corr']
df = df.sort_values(by = ['loc_ssize', 'word_ssize'])

df.to_csv('gtrp_diareg_overlap_simulation.txt', sep = '\t', index=False)
