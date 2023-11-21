import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def plot_inferred(df, method):
    probas = np.array(sorted(set(df["proba_df"])))

    prob_inf = []
    prob_std = []
    for proba in probas:
        prob_inf += [df[df["proba_df"] == proba]["probas_inferred"].mean()]
        prob_std += [df[df["proba_df"] == proba]["probas_inferred"].std()]

    prob_inf = np.array(prob_inf)
    prob_std = np.array(prob_std)

    if method in {"beast", "real"}:
        rate_inf = []
        rate_std = []
        for proba in probas:
            rate_inf += [df[df["proba_df"] == proba]["rates_inferred"].mean()]
            rate_std += [df[df["proba_df"] == proba]["rates_inferred"].std()]
            
        rate_inf = np.array(rate_inf)
        rate_std = np.array(rate_std)

        fig, axs = plt.subplots(2, 1)

        axs[0].plot(probas, prob_inf)
        axs[0].plot(probas, probas)
        axs[0].fill_between(probas, prob_inf - 1.96*prob_std, prob_inf + 1.96*prob_std, color='b', alpha=.1)
        axs[0].plot(probas, 0.5*(1 - np.exp(-2*rate_inf)))

        axs[1].plot(probas, rate_inf)
        axs[1].plot(probas, -np.log(1-2*probas)/2)
        axs[1].fill_between(probas, rate_inf - 1.96*rate_std, rate_inf + 1.96*rate_std, color='b', alpha=.1)

        axs[0].set_xlabel('$p_{true}$')
        axs[0].title.set_text("$p$ inferred")
        axs[1].set_xlabel('$p_{true}$')
        axs[1].title.set_text("rate inferred")

        axs[0].legend(['$p_{inferred}$', '$p_{true}$', '0.95 confidene interval', '$\\frac{1}{2}(1 - \exp(-2 \lambda_{inferred}))$'])
        axs[1].legend(['$\lambda_{inferred}$', '$-\\frac{1}{2}log(1-2p)$'])

        fig.suptitle(f'Number of different trees is {len(set(df["exp_df"]))}, length of MCMC chain is 10000, method is {method}', fontsize=16)
        fig.set_figheight(12)
        fig.set_figwidth(15)
    else:
        fig, axs = plt.subplots(1,1)

        axs.plot(probas, prob_inf)
        axs.plot(probas, probas)
        axs.fill_between(probas, prob_inf - 1.96*prob_std, prob_inf + 1.96*prob_std, color='b', alpha=.1)

        axs.set_xlabel('$p_{true}$')
        axs.title.set_text("$p$ inferred")

        axs.legend(['$p_{inferred}$', '$p_{true}$', '0.95 confidene interval'])

        fig.suptitle(f'Number of different trees is {len(set(df["exp_df"]))}, method is {method}', fontsize=16)
        fig.set_figheight(12)
        fig.set_figwidth(15)

    fig.savefig(f'plot_experiment_{method}.jpg')

def plot_inferred_weighted(df, method):
    probas = np.array(sorted(set(df["proba_df"])))
    
    prob_inf = []
    prob_01_inf = []
    prob_10_inf = []
    for proba in probas:
        prob_inf += [df[df["proba_df"] == proba]["probas_inferred"].mean()]
        prob_01_inf += [df[df["proba_df"] == proba]["weighted_proba_01"].mean()]
        prob_10_inf += [df[df["proba_df"] == proba]["weighted_proba_10"].mean()]
        
    prob_inf = np.array(prob_inf)
    prob_01_inf = np.array(prob_01_inf)
    prob_10_inf = np.array(prob_10_inf)

    fig, axs = plt.subplots(1,1)

    axs.plot(probas, prob_inf)
    axs.plot(probas, prob_01_inf)
    axs.plot(probas, prob_10_inf)
    axs.plot(probas, probas)

    axs.set_xlabel('$p_{true}$')
    axs.title.set_text("$p$ inferred")

    axs.legend(['$p_{inferred}$', '$p_{01}$ inferred', '$p_{10}$ inferred',  '$p_{true}$'])

    fig.suptitle(f'Number of different trees is {len(set(df["exp_df"]))}, method is {method}', fontsize=16)
    fig.set_figheight(12)
    fig.set_figwidth(15)

    fig.savefig(f'plot_experiment_weighted_{method}.jpg')

def plot_one_figure(dfs, methods):
    same_probas = True
    probas = set(dfs[0]["proba_df"])
    for i, df in enumerate(dfs):
        probas_new = set(df["proba_df"])
        if probas != probas_new:
            same_probas = False
            break
        probas = probas_new
    #if same_probas
    probas = np.array(list(sorted(probas)))
    fig, axs = plt.subplots(1,1)

    axs.plot(probas, probas)
    for i, df in enumerate(dfs):
        prob_inf = []
        for proba in probas:
            prob_inf += [df[df["proba_df"] == proba]["probas_inferred"].mean()]

        prob_inf = np.array(prob_inf)

        axs.plot(probas, prob_inf)
    
    axs.set_xlabel('$p_{true}$')
    axs.title.set_text("$p$ inferred")
    axs.legend(["ground truth"] + methods)
    fig.suptitle(f'Number of different trees is {len(set(df["exp_df"]))}', fontsize=16)
    fig.set_figheight(12)
    fig.set_figwidth(15)
    fig.savefig(f'plot_all_methods_proba.jpg')

    methods_rate = []
    cols = [set(df.columns) for df in dfs]
    col_set = set.union(*cols)
    if "rates_inferred" in col_set:
        fig, axs = plt.subplots(1,1)
        
        axs.plot(probas, -np.log(1-2*probas)/2)
        for i, df in enumerate(dfs):
            if "rates_inferred" in cols[i]:
                methods_rate += methods[i],
                rate_inf = []
                for proba in probas:
                    rate_inf += [df[df["proba_df"] == proba]["rates_inferred"].mean()]

                rate_inf = np.array(rate_inf)
                axs.plot(probas, rate_inf)

        axs.set_xlabel('$p_{true}$')
        axs.title.set_text("rate inferred")
        axs.legend(["ground truth"] + methods_rate)
        fig.suptitle(f'Number of different trees is {len(set(df["exp_df"]))}', fontsize=16)
        fig.set_figheight(12)
        fig.set_figwidth(15)
        fig.savefig(f'plot_all_methods_rates.jpg')

def plot_errors(dfs, methods):
    same_probas = True
    probas = set(dfs[0]["proba_df"])
    for i, df in enumerate(dfs):
        probas_new = set(df["proba_df"])
        if probas != probas_new:
            same_probas = False
            break
        probas = probas_new
    #if same_probas
    probas = np.array(list(sorted(probas)))


    cols = [set(df.columns) for df in dfs]
    col_set = set.union(*cols)

    fig_with_rate, axs_with_rate = plt.subplots(1,1)
    fig_no_rate, axs_no_rate = plt.subplots(1,1)

    methods_no_rate = []
    methods_with_rate = []
    for i, df in enumerate(dfs):
        prob_inf = []
        for proba in probas:
            prob_inf += [df[df["proba_df"] == proba]["probas_inferred"].mean()]

        prob_inf = np.array(prob_inf)

        if "rates_inferred" in df.columns:
            methods_with_rate += methods[i],
            axs_with_rate.plot(probas, prob_inf - probas)
        else:
            methods_no_rate += methods[i],
            axs_no_rate.plot(probas, prob_inf - probas)

    axs_no_rate.plot(probas, np.zeros(len(probas)), "r", linewidth = 2)
    axs_no_rate.set_xlabel('$p_{true}$')
    axs_no_rate.title.set_text("error $p$ inferred")
    axs_no_rate.legend(methods_no_rate + ["0"])
    fig_no_rate.suptitle(f'Number of different trees is {len(set(df["exp_df"]))}', fontsize=16)
    fig_no_rate.set_figheight(12)
    fig_no_rate.set_figwidth(15)
    fig_no_rate.savefig(f'plot_proba_errors_without_rate.jpg')

    axs_with_rate.plot(probas, np.zeros(len(probas)), "r", linewidth = 2)
    axs_with_rate.set_xlabel('$p_{true}$')
    axs_with_rate.title.set_text("error $p$ inferred")
    axs_with_rate.legend(methods_with_rate + ["0"])
    fig_with_rate.suptitle(f'Number of different trees is {len(set(df["exp_df"]))}', fontsize=16)
    fig_with_rate.set_figheight(12)
    fig_with_rate.set_figwidth(15)
    fig_with_rate.savefig(f'plot_proba_errors_with_rate.jpg')
    
    methods_rate = []
    cols = [set(df.columns) for df in dfs]
    col_set = set.union(*cols)
    if "rates_inferred" in col_set:
        fig, axs = plt.subplots(1,1)
        for i, df in enumerate(dfs):
            if "rates_inferred" in df.columns:
                methods_rate += methods[i],
                rate_inf = []
                for proba in probas:
                    rate_inf += [df[df["proba_df"] == proba]["rates_inferred"].mean()]

                rate_inf = np.array(rate_inf)

                axs.plot(probas, rate_inf + np.log(1-2*probas)/2)
        axs.plot(probas, np.zeros(len(probas)), "r", linewidth = 2)
        axs.set_xlabel('$p_{true}$')
        axs.title.set_text("error rate inferred")
        axs.legend(methods_rate + ["0"])
        fig.suptitle(f'Number of different trees is {len(set(df["exp_df"]))}', fontsize=16)
        fig.set_figheight(12)
        fig.set_figwidth(15)
        fig.savefig(f'plot_rate_errors_with_rate.jpg')


df_real = pd.read_csv("experiment_real/output_real_tree_accurracy.csv")
df_beast = pd.read_csv("experiment_beast/output_beast_tree_experiment.csv")
df_heur = pd.read_csv("experiment_heuristic/df_output_heuristic_experiment.csv")
df_heur_informed = pd.read_csv("experiment_heuristic_informed/df_output_heuristic_informed_experiment.csv")

dfs = [df_real, df_beast, df_heur, df_heur_informed]
methods = ["baseline_real_tree", "baseline_beast_tree", "heuristic", "heuristic_informed"]

plot_inferred(df_beast, method = "baseline_beast_tree")
plot_inferred(df_real, method = "baseline_real_tree")
plot_inferred(df_heur, method = "heuristic")
plot_inferred(df_heur_informed, method = "heuristic_informed")

plot_errors(dfs, methods)

plot_one_figure(dfs, methods)

plot_inferred_weighted(df_beast, method = "baseline_beast_tree")
plot_inferred_weighted(df_real, method = "baseline_real_tree")