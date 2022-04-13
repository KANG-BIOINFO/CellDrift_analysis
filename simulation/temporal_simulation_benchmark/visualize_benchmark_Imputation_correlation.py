import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib as mpl

# global parameters
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

def create_cmap(cmap_name, cmap_list):
    cmap_o = plt.get_cmap(cmap_name)
    for i in range(cmap_o.N):
        rgba = cmap_o(i)
        cmap_list.append(mpl.colors.rgb2hex(rgba))
    return cmap_list

cmap = []; cmap = create_cmap('tab10', cmap)

### 
df_results_kmeans = pd.read_csv('results_clusters_kmeans_fuzzy.txt', sep = '\t', header = 0)
df_results_kmeans_DEGs = pd.read_csv('results_clusters_kmeans_fuzzy_DEGs.txt', sep = '\t', header = 0)
df_results_EMCluster = pd.read_csv('results_clusters_fPCA_EMCluster.txt', sep = '\t', header = 0)
df_results_DEGs_EMCluster = pd.read_csv('results_clusters_DEG_fPCA_EMCluster.txt', sep = '\t', header = 0)

df_metrics = pd.concat([df_results_kmeans, df_results_kmeans_DEGs, df_results_EMCluster, df_results_DEGs_EMCluster], axis = 0)

covariate_names = df_metrics['covariate_name'].unique()

# n_rep = 3
n_rep = 12
n_methods = 3
colors = [cmap[3]] + cmap[5:7]
alphas = [1, 0.15, 0.15]
methods = ['FPCA_EMCluster', 'FPCA_EMCluster_t-test', 'FPCA_EMCluster_wilcoxon']

for covariate_name in covariate_names:
    df_sub = df_metrics.loc[df_metrics['covariate_name'] == covariate_name, : ].copy()
    df_mean = df_sub.groupby(['method', 'covariate_level']).mean()['correlation']
    df_std = df_sub.groupby(['method', 'covariate_level']).std()['correlation'] / np.sqrt(n_rep)

    fig, ax = plt.subplots()
    for method, color, alpha in zip(methods, colors, alphas):
        # get mean
        df_mean2 = df_mean.loc[method, :].copy()

        covariate_values = [i[1] for i in df_mean2.index.values]
        correlation_score_means = df_mean2.values
        ax.plot(covariate_values, correlation_score_means, label = method, color = color, alpha = alpha)
        ax.scatter(covariate_values, correlation_score_means, color = color, alpha = alpha)

        # get standard errors
        df_std2 = df_std.loc[method, :].copy()
        covariate_values = [i[1] for i in df_std2.index.values]
        correlation_score_errors = df_std2.values
        ax.errorbar(x = covariate_values, y = correlation_score_means, yerr = correlation_score_errors, color = color, alpha = alpha)
        ax.grid(True, alpha = 0.2)

    ax.set_ylim(0.7, 1.05)
    handles, labels = ax.get_legend_handles_labels()
    plt.legend(handles = handles, labels = labels)
    fig.savefig('figures/benchmark_imputation_' + covariate_name + '.pdf', bbox_inches = 'tight')

