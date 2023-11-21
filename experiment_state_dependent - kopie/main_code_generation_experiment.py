# %%
import multiprocessing
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy.sparse.linalg as lin
import scipy.sparse as sp
import pandas as pd
import networkx as nx
from networkx.readwrite import json_graph
import time
import json
import os
import copy
import functools
import itertools
import pickle
import random

from dotdict import dotdict

from agent_model import Agent_Country


def num2gen(arr):
    s=""
    for i in arr:
        if i==0:
            s+="G"
        elif i==1:
            s+="A"
        elif i==2:
            s+="C"
        elif i==3:
            s+="T"
    return s

def num2state(num):
    if num==0:
        return "A"
    elif num==1:
        return "B"
    else:
        return "C"

def run_one_epi(graph_name, beta, n, it, distribution, probas, gamma=1., n_sites = 10, sub_rate = 0.5):
    np.random.seed(it)

    if graph_name == "complete":
        graph_args = dotdict({
            "name":"complete",
            "N": n,
        })
    else:
        raise Exception("graph_name not understood")


    #print(graph_args)
    graph=init_graph(graph_args)

    args = dotdict({
        "plot": False,
        "max_iteration": 10000,
        "beta": beta,
        "beta_super":0.0,
        "sigma": 0.5,
        "gamma": gamma,
        "xi":0,
        "p_teleport":0.0,
        "MAX_E_TIME":10,
        "MAX_I_TIME":10,
        "infected_agents": [0],
        "super_infected_agents": [],
        "p_super": 0.0,
        "distribution": distribution,
        "probas": probas
    })
    args["I_time"]=int(1/args["gamma"])
    #print(args)

    country = Agent_Country(args, graph)
    num_inf=[]
    #country.log_json()
    for i in range(args["max_iteration"]):
        #results[-1].append(country.states==2)
        #temp_results.append(np.sum(country.states==2))
        if country.check_stop():
            #print("Infection extinction at iteration: {}".format(i))
            break
        num_inf.append(np.sum(country.states==2))
        country.step()
        #country.log_json()
    #results[-1].append(country.states==3)

    names = {}
    parent = {}
    tree = nx.DiGraph()
    for name in country.graph.nodes():
        if "infector" in country.graph.nodes[name]:
            if country.graph.nodes[name]["infector"]>=0:
                tree.add_edge(country.graph.nodes[name]["infector"], country.graph.nodes[name]["index"])
                parent[country.graph.nodes[name]["index"]] = country.graph.nodes[name]["infector"]
            names[graph.nodes[name]["index"]]=name

    nx.set_node_attributes(tree, parent, 'parent')
    codes_dict={}
    times_dict={}
    for i in nx.topological_sort(tree):
        times_dict[i]=country.graph.nodes[names[i]]["first_inf_time"]

        infector = country.graph.nodes[names[i]]["infector"]
        if infector==-1:
            codes_dict[i]=np.zeros(n_sites)
        else:
            # parent gcode
            codes_dict[i]=copy.deepcopy(codes_dict[infector])
            # choosing substitution indices
            sub_ind = np.random.choice(np.arange(n_sites), np.random.binomial(n_sites, sub_rate/n_sites ))
            # making the substitution
            codes_dict[i][sub_ind]=np.mod(codes_dict[i][sub_ind]+np.random.randint(1,4,size=len(sub_ind)),4)
    
    time_df=pd.DataFrame(np.cumsum(num_inf)).transpose()
    time_df.columns=np.array(time_df.columns)
    
    df=pd.DataFrame.from_dict(nx.get_node_attributes(country.graph,"first_inf_time"),orient='index')
    df.columns=["first_inf_time"]
    df["index"]=df.index
    df=df[["index","first_inf_time"]]
    df=pd.merge(df,pd.DataFrame(pd.Series(codes_dict).map(num2gen)), left_index=True, right_index=True)
    df.columns=["strain","date","code"]
    states = [country.type_dict[key] for key in df["strain"]]
    df["state"] = states

    return country, tree, names, df, time_df

def init_graph(args):
    if (args.name == "complete"):
        graph=nx.complete_graph(args.N)
        nx.set_node_attributes(graph, {n:i for i,n in enumerate(graph.nodes())}, "index")
    else:
        print("Graph mode not understood", args.graph)

    return graph

def write_to_fasta(df, seed):
    with open(f"experiment_beast/code_exp_{seed}.fasta", "w") as file:
        for row in df.iterrows():
            file.write(">"+str(row[1]["strain"])+" "+str(row[1]["date"])+" "+str(row[1]["state"])+"\n")
            file.write(row[1]["code"]+"\n")

def reconstruct_tree(tips, tree):
    all_paths = [nx.shortest_path(tree, 0, leaf) for leaf in tips]
    nodes = set(itertools.chain.from_iterable(all_paths))
    return nx.subgraph(tree, nodes)

def DiTree_to_newick_string(subtree, node):
    childeren = list(subtree.neighbors(node))
    if len(childeren) == 0:
        return str(node)+":1"
    else:
        return "(" + "".join(str(DiTree_to_newick_string(subtree, child)) + \
                             "," for child in childeren[:-1])+ str(DiTree_to_newick_string(subtree, childeren[-1])) + "):1"

def DiTree_to_newick_string_bifur(subtree, node, correction):
    childeren = list(subtree.neighbors(node))
    if len(childeren) == 0:
        return str(node)+f":{1-correction}"
    elif len(childeren) == 1:
        return str(DiTree_to_newick_string_bifur(subtree, childeren[0], correction - 1))

    elif len(childeren) == 2:
        return "(" + str(DiTree_to_newick_string_bifur(subtree, childeren[0], 0)) + "," + str(DiTree_to_newick_string_bifur(subtree, childeren[1], 0)) + f"):{1-correction}"
    else:
        random.shuffle(childeren)
        corrections = [random.random()/100 for _ in range(len(childeren) - 2)]
        string = "(" + DiTree_to_newick_string_bifur(subtree, childeren.pop(), sum(corrections)) + "," + DiTree_to_newick_string_bifur(subtree, childeren.pop(), sum(corrections)) + f"):{corrections[-1]}"
        for i, child in enumerate(childeren[:-1]):
            string = "(" + DiTree_to_newick_string_bifur(subtree, child, sum(corrections[:-(1+i)])) + "," + string + f"):{corrections[-(2+i)]}"
        return "(" + DiTree_to_newick_string_bifur(subtree, childeren[-1], 0) + "," + string + f"):{1-correction}"

def write_to_nexus(newick_string, ntax, df, seed, matrix_i):
    with open(f"experiment_real/real_tree_{matrix_i}_exp_{seed}.nex", "w") as file:
        file.write(f"#NEXUS\nBegin TAXA;\n  Dimensions ntax={ntax};\n  TaxLabels\n")
        for  _, row in df.iterrows():
            file.write("  '"+str(row["strain"])+" "+str(row["date"])+" "+str(row["state"])+"'\n")
        file.write("  ;\n")
        file.write("End;\n\nBegin trees;\n  Translate\n")
        counter = 0
        for _, row in df.iterrows():
            counter += 1
            file.write(f"  {counter} '"+str(row["strain"])+" "+str(row["date"])+" "+str(row["state"])+"',\n")
        file.write("  ;\n")
        file.write("  tree TREE1 = [&R] " + newick_string + ";\nEND;")

def run_and_write_nexus(df, tree, seed, matrix_i, fasta = False, nexus = False):
    if fasta:
        write_to_fasta(df, seed)
    if nexus:
        tips = set(df["strain"])
        ntax = len(tips)
        subtree = reconstruct_tree(tips, tree) # also this has to happen before the loop inside an "if nexus" statement

        # relabel the nodes such that all the leaves are labeled 0 to ntax-1
        # otherwise the ape (R library) function raises an error
        mapping_1 = {strain: strain_i + 1 for strain_i, strain in enumerate(df["strain"])}
        previous_nodes = set(subtree.nodes()).difference(tips)
        mapping_2 = {strain: strain_i + ntax + 1 for strain_i, strain in enumerate(sorted(previous_nodes))}
        mapping = {**mapping_1, **mapping_2}

        subtree_relabled = nx.relabel_nodes(subtree, mapping)

        newick_string = DiTree_to_newick_string_bifur(subtree_relabled, ntax + 1, correction = 0)
        # up to this can go before the loop
        write_to_nexus(newick_string, ntax, df, seed, matrix_i)

def write_df_childeren_to_csv(tree, df, seed, heuristic=False, heuristic_informed=False):
    if heuristic:
        df.to_csv(f"experiment_heuristic/df_childeren_exp_{seed}.csv", index = False)
    if heuristic_informed:
        tree_df = nx.to_pandas_edgelist(tree)
        tree_df.columns = ["parent", "child"]
        tree_df = tree_df.iloc[:, ::-1]
        tree_df = tree_df.sort_values(by=['child'])
        tree_df.to_csv(f"experiment_heuristic_informed/df_tree_exp_{seed}.csv", index = False)

def write_df_types_to_csv(df, proba, seed):
    df.to_csv(f"experiment_heuristic/df_types_{int(proba*100)}_exp_{seed}.csv", index = False)



def main(graph_name = "complete", n = 500, distributions = None, contact_matrices = None, beta = 0.5, num_seeds = 100, sampling = 1, gamma = 1, fasta=False, nexus=False, heuristic=False, heuristic_informed=False):
    # to do:
    # write_to_fasta needs to write a .nex file instead of .fasta so we don't have to use beastgen with fasta_to_nexus.template
    # rename write_to_fasta
    # use \t instead of spaces
    seed_and_matrix_i = []
    for matrix_i, contact_matrix in enumerate(contact_matrices):
        distribution = distributions[matrix_i]
        seeds = num_seeds*[None]
        if distribution == None:
            distribution = [int(n/2), n - int(n/2)]
            row_product = [1/k if k > 0 else 0 for k in distribution] # 0 individuals in group k means zero probability
            transition_matrix = beta*np.array([[5, 5], [5, 5]]).dot(np.diag(row_product))
        elif sum(distribution) != n:
            print("Population distribution does sum not add up to n")
        elif max(distribution) == n:
            print("Population can't be contained in one group")
        elif contact_matrix.max() == 0:
            print("Contact matrix can't be zero")
        elif contact_matrix.shape != (len(distribution), len(distribution)):
            print("Dimensions of the contact-matrix does not match the number of groups")
        elif np.amax(contact_matrix)*beta/n > 1:
            print("Resulting P-matrix will contain an element with p_ij > 1")
        else:
            row_product = [1/k if k > 0 else 0 for k in distribution] # 0 individuals in group k means zero probability
            transition_matrix = beta*contact_matrix.dot(np.diag(row_product))
            seed, seeds = 0, []
        
        while len(seeds) < num_seeds:
            seed += 1
            _, tree, _, df, _ = run_one_epi(graph_name, beta, n, seed, distribution, transition_matrix, gamma = gamma)
            df_sample = df[df["date"] == df["date"].mode().max()]
            state_set = set(df_sample["state"])
            if len(df_sample) > 1 and len(state_set) > 1:
                seeds += seed,
                seed_and_matrix_i += [[seed, matrix_i]]
                df_sample = df[df["date"] == df["date"].mode().max()]
                run_and_write_nexus(df_sample, tree, seed, matrix_i, fasta=fasta, nexus=nexus)


    if distributions == None:
        with open("contact_matrix.csv", "w") as file:
            file.write(f"None")
        with open("distribution.csv", "w") as file:
            file.write(f"None")
    else:
        np.savetxt("contact_matrix.csv", contact_matrix, delimiter=",")
        np.savetxt("distribution.csv", np.array(distribution), delimiter=",")
        with open("matrices.txt", "w") as file:
            for matrix_i, contact_matrix in enumerate(contact_matrices):
                file.write(str(matrix_i) + " " + str(contact_matrix.tolist()) + "\n")
        with open("distributions.txt", "w") as file:
            for distribution_i, distribution in enumerate(distributions):
                file.write(str(distribution_i) + " " + str(distribution) + "\n")

    with open("experiments_and_matrix_indices.txt", "w") as file:
        for experiment in seed_and_matrix_i:
            file.write(str(experiment[0]) + " " + str(experiment[1]) + "\n")

    with open("parameters.txt", "w") as file:
        file.write(f"beta,gamma,graph_name,graph_size,sampling\n")
        file.write(f"{beta},{gamma},{graph_name},{n},{sampling}")


if __name__ == "__main__":
    contact_matrices = [np.array([[x, x], [x, x]]) for x in range(17, 20, 2)]
    distributions = len(contact_matrices)*[[250, 250]]
    main(n = 500, distributions = distributions, contact_matrices = contact_matrices, num_seeds=4, nexus = True)