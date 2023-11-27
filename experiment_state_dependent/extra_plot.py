import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import itertools


df_real = pd.read_csv("experiment_real/output_real_tree_BiSSE.csv")

exps_file = open('experiments.txt', 'r')
exps_lines = exps_file.readlines()
exps = [int(exp_line[:-1]) for exp_line in exps_lines]


rate_00 = df_real["rate_00"].mean()
rate_01 = df_real["rate_01"].mean()
rate_10 = df_real["rate_10"].mean()
rate_11 = df_real["rate_11"].mean()

for rate in df_real.columns[1:]:
    print(f"{rate} = {0.5*(1-np.exp(-2*np.array(df_real[rate].mean())))}")