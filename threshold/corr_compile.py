import natsort
import pandas as pd
import pathlib

cor_vector_list = [path for path in pathlib.Path('./matrices').glob('*.txt')]

list_of_cor_vecs = {}
for file in natsort.natsorted(cor_vector_list, alg=natsort.ns.PATH):
	with open(file, 'r') as read_file:
		# print()
		# print(read_file.readline())
		list_of_cor_vecs.update({str(file.name).strip('.txt').strip('aggregate_matrix_'): read_file.readline().split('\t')})

# print(list_of_cor_vecs)

df = pd.DataFrame.from_dict(list_of_cor_vecs)
print(df)
# df.columns = ['ss_' + str(i) for i in range(1,101)]

df.to_csv('cor_vec_gtrp.txt', sep = '\t', index=False)