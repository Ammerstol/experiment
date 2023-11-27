import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

country = "belgium"
week = "2"
df = pd.read_csv(f"gisaid_data/{country}_wk_{week}_metadata.tsv", sep="\t")
if df["age"].dtype == "object":
    ages_list = [str(i) for i in range(117)] + ["1 month"] + [f"{i} months" for i in range(2, 11)]
    df.drop(df[~df["age"].isin(ages_list)].index, inplace = True)
    df["age_category"] = list(map(len, df["age"]))

    df.loc[df["age_category"] > 3, "age"] = 0
    df["age"] = df["age"].astype(int)
elif df["age"].dtype == "float64" or df["age"].dtype == "int64":
    df = df[df["age"].notnull()]
    df["age"] = df["age"].astype(int)
else:
    print("data in incorrect format")

print(len(df))
df.drop_duplicates(subset = "strain", keep = 'first', inplace = True)
print(len(df))

multi_fasta_file = open(f"gisaid_data/{country}_wk_{week}_sequences.fasta")
multi_lines = multi_fasta_file.readlines()

if True:
    new_lines = []
    strain_set = set(df["strain"])
    write_to_file = False
    for line in multi_lines:
        if line[0] == ">":
            if line[1:-1] in strain_set:
                write_to_file = True
            else:
                write_to_file = False
        if write_to_file:
            new_lines += line,

    new_fasta_file = open(f"gisaid_data/{country}_wk_{week}_filtered_sequences.fasta", 'w')
    new_fasta_file.writelines(new_lines)
    new_fasta_file.close()

if False:
    df2 = df.pivot_table(index = ['strain'], aggfunc ='size')
    print(list(df2).count(2))

if False:
    plt.hist(df["age"], bins = range(0, max(df["age"]) + 10, 10))
    plt.show()
