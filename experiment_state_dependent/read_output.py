import pandas as pd
import numpy as np
from scipy import linalg
import csv
import os
import matplotlib.pyplot as plt

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

params = pd.read_csv("parameters.txt")
beta = params["beta"].item()

df_MuSSE = pd.read_csv("experiment_real/output_real_tree_MuSSE.csv")

cols = df_MuSSE.columns
transition_rates = sorted([name for name in cols if name[0] == "q"])
speciation_rates = sorted([name for name in cols if name[0] == "l"])
all_rates = speciation_rates + transition_rates

val_dict = {exp: {matrix_i: [] for matrix_i in matrix_indices} for exp in all_rates}

for row_i, row in df_MuSSE.iterrows():
    for exp in all_rates:
        val_dict[exp][mats_exps[row_i][0]] += row[exp],

if True:
    df_bd = pd.read_csv("experiment_real/output_real_tree_bd.csv")

    cols = df_bd.columns
    speciation_rates = sorted([name for name in cols if name[0] == "l"])

    val_dict_bd = {matrix_i: [] for matrix_i in matrix_indices}

    for row_i, row in df_bd.iterrows():
            val_dict_bd[mats_exps[row_i][0]] += row["lambda"],

print(val_dict_bd)

k = int(np.ceil(np.sqrt(len(all_rates))))
fig, axes = plt.subplots(k, k)

for ax_i, ax in enumerate(axes.flat):
    rate_list = []
    rate_list_bd = []
    for matrix_i in matrix_indices:
        rate_list += np.array(val_dict[all_rates[ax_i]][matrix_i]).mean(),
        if all_rates[ax_i][0] == "l":
            rate_list_bd += np.array(val_dict_bd[matrix_i]).mean(),
    ax.plot(rate_list, label = "MuSSE")
    if all_rates[ax_i][0] == "l":
        ax.plot(rate_list_bd, label = "BD")
        handles, labels = ax.get_legend_handles_labels()
    ax.set_title(all_rates[ax_i])
fig.legend(handles, labels, loc=(0, 0.84))
fig.suptitle("Comparison between MuSSE and Birth-Death (BD) parameter estimation,\nnote there is only one lambda in the BD process")
fig.tight_layout()
plt.show()


















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