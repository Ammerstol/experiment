import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import itertools
import os

os.chdir(os.path.dirname(os.path.abspath(__file__)))

df_bd = pd.read_csv("experiment_real/output_real_tree_bd_single.csv")
plt.hist(df_bd["lambda"], bins = 40)
plt.show()