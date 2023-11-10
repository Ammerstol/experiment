import subprocess
import fasta_plus_mcc_to_mcc
import main_code_generation_experiment
import pandas as pd
import numpy as np
import time
import itertools

def estimate_parents(df, gamma):
    parents = []
    max_lik = []
    for child in df.iterrows():
        w = [compute_weight(child, row, gamma) for row in df.iterrows()]
        max_lik += np.amax(w),
        arg_max_array = np.argwhere(w == max_lik[-1]).flatten()
        parents += np.random.choice(arg_max_array),
    # also return a dict of weights for a weighted p estimation?
    return parents, max_lik

def compute_weight(child, row, gamma):
    time_diff = child[1]["date"] - row[1]["date"]
    base_diff = compare_sites(child[1]["code"], row[1]["code"])
    if time_diff <= 0:
        return 0
    else:
        return gamma**time_diff * 0.25**base_diff

def compare_sites(child_code, row_code):
    base_diff = 0
    for child_base, row_base in zip(child_code, row_code):
        if child_base != row_base:
            base_diff += 1

    return base_diff

def estimate_jump_proba(parent_types, childeren_types):
    return np.mean(np.abs(parent_types - childeren_types))


heuristic = True
heuristic_informed = True

start = time.time()
main_code_generation_experiment.main(num_seeds=100, heuristic=heuristic, heuristic_informed=heuristic_informed)


probas_file = open('true_probabilities.txt', 'r')
probas_lines = probas_file.readlines()
probas = [int(float(proba_line[:-1])*100) for proba_line in probas_lines]

exps_file = open('experiments.txt', 'r')
exps_lines = exps_file.readlines()
exps = [int(exp_line[:-1]) for exp_line in exps_lines]

params_df = pd.read_csv("parameters.txt")

loop = time.time()
print(f"code generation is {loop - start}")

proba_infer_heur = []
proba_infer_heur_informed = []

for exp in exps:
    df_sample_hear = pd.read_csv(f"experiment_heuristic/df_childeren_exp_{exp}.csv")
    if heuristic:
        df_sample_hear = pd.read_csv(f"experiment_heuristic/df_childeren_exp_{exp}.csv")
        np.random.seed(exp)
        parents_index_heur, _ = estimate_parents(df_sample_hear, params_df["gamma"].item())
    if heuristic_informed:
        df_real_tree = pd.read_csv(f"experiment_heuristic_informed/df_tree_exp_{exp}.csv")
        list_of_index_sets = list(df_real_tree[df_real_tree["child"] == parent].index for parent in df_real_tree["parent"])
        parents_index_heur_informed = [index_set[0]+1 if len(index_set)>0 else 0 for index_set in list_of_index_sets]
    for proba in probas:
        childeren_types = np.array(pd.read_csv(f"experiment_heuristic/df_types_{proba}_exp_{exp}.csv")["state"])
        if heuristic:
            parent_types_heur = np.array([childeren_types[i] for i in parents_index_heur])
            proba_infer_heur += estimate_jump_proba(parent_types_heur[1:], childeren_types[1:]),
        if heuristic_informed:
            parent_types_heur_informed = np.array([childeren_types[i] for i in parents_index_heur_informed])
            proba_infer_heur_informed += estimate_jump_proba(parent_types_heur_informed, childeren_types[1:]),

        
end = time.time()
print(f"looping is {end - loop}")
exp_df = list(itertools.chain.from_iterable(len(probas)*[exp] for exp in exps))
proba_df = list(itertools.chain.from_iterable(np.array(probas)/100 for _ in exps))
if heuristic:
    df_output = pd.DataFrame(data={"exp_df": exp_df,
                            "proba_df": proba_df,
                            "probas_inferred": proba_infer_heur})
    df_output.to_csv(f"experiment_heuristic/df_output_heuristic_experiment.csv", index = False)
if heuristic_informed:
    df_output = pd.DataFrame(data={"exp_df": exp_df,
                            "proba_df": proba_df,
                            "probas_inferred": proba_infer_heur_informed})
    df_output.to_csv(f"experiment_heuristic_informed/df_output_heuristic_informed_experiment.csv", index = False)