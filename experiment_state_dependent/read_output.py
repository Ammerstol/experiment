import pandas as pd
import numpy as np
from scipy import linalg
import csv
import os

os.chdir(os.path.dirname(os.path.abspath(__file__)))

with open('matrices.txt', 'r') as f:
    csv_file = csv.reader(f,delimiter=' ')
    matrix_indices = [int(row[0]) for row in csv_file]
    matrices_delimited = [row[1:] for row in csv_file]

with open('experiments_and_matrix_indices.txt', 'r') as f:
    lines = f.readlines()
    mats_exps_str = [line.split(" ") for line in lines]
    mats_exps = [(int(str_list[0]), int(str_list[1])) for str_list in mats_exps_str]

matrices = []
for delimited in matrices_delimited:
    matrices += [np.array(eval(''.join(delimited)))]

val_dict = {matrix_i: [] for matrix_i in matrix_indices}

params = pd.read_csv("parameters.txt")
beta = params["beta"].item()

df_real = pd.read_csv("experiment_real/output_real_tree_bd.csv")

for row_i, row in df_real.iterrows():
    val_dict[mats_exps[row_i][0]] += row["lambda"],

for key in val_dict:
    print(f"{key}: {sum(val_dict[key])/len(val_dict[key])}")

# for mat_exp in 

# cols = df_real.columns
# transition_rates = sorted([name for name in cols if name[0] == "q"])
# speciation_rates = sorted([name for name in cols if name[0] == "l"])
# dim = max([int(rate[-1]) for rate in transition_rates])

# matrix_indices = list(sorted(set(df_real["matrix_i"])))
# for matrix_i in matrix_indices:
#     df_sample = df_real[df_real["matrix_i"] == matrix_i]

#     speciation_vector = np.zeros((dim))
#     intensity_matrix = np.zeros((dim, dim))

#     for rate in transition_rates:
#         intensity_matrix[int(rate[-2])-1, int(rate[-1]) - 1] = df_sample[rate].mean()

#     intensity_matrix = intensity_matrix - np.diag([sum(rate_row) for rate_row in intensity_matrix])
#     transition_matrix_1 = linalg.expm(intensity_matrix)

#     for rate in speciation_rates:
#         speciation_vector[int(rate[-1])-1] = df_sample[rate].mean()
    
#     speciation_vector = speciation_vector # np.exp(speciation_vector)
#     contact_matrix_inferred = np.diag(speciation_vector).dot(transition_matrix_1)

#     beta_x_contact_matrix = beta*matrices[matrix_i]
#     error = linalg.norm(beta_x_contact_matrix - contact_matrix_inferred)

#     print(f"Frobenius error for index {matrix_i}:\n{error}\n")
#     print(speciation_vector)
#     print(intensity_matrix)
#     print(transition_matrix_1)
#     print(contact_matrix_inferred)
#     print(beta_x_contact_matrix)
#     print("")