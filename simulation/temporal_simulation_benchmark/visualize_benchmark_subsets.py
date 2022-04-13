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

n_rep = 3
n_methods = 5
colors = cmap[: n_methods]
methods = ['kmeans', 'kmeans_fuzzy_DEGs_t-test', 'kmeans_fuzzy_DEGs_wilcoxon', 'kmeans_fuzzy', 'FPCA_EMCluster']
alphas = [0.2, 0.2, 0.2, 1, 0.2]

for covariate_name in covariate_names:
    df_sub = df_metrics.loc[df_metrics['covariate_name'] == covariate_name, : ].copy()
    df_mean = df_sub.groupby(['method', 'covariate_level']).mean()['ARI']
    df_std = df_sub.groupby(['method', 'covariate_level']).std()['ARI'] / np.sqrt(n_rep)

    fig, ax = plt.subplots()
    for method, color, alpha in zip(methods, colors, alphas):
        # get mean
        df_mean2 = df_mean.loc[method, :].copy()

        covariate_values = [i[1] for i in df_mean2.index.values]
        ARI_score_means = df_mean2.values
        ax.plot(covariate_values, ARI_score_means, label = method, color = color, alpha = alpha)
        ax.scatter(covariate_values, ARI_score_means, color = color, alpha = alpha )

        # get standard errors
        df_std2 = df_std.loc[method, :].copy()
        covariate_values = [i[1] for i in df_std2.index.values]
        ARI_score_errors = df_std2.values
        ax.errorbar(x = covariate_values, y = ARI_score_means, yerr = ARI_score_errors, color = color, alpha = alpha)
        ax.grid(True, alpha = 0.2)

    ax.set_ylim(0.2, 1.1)
    handles, labels = ax.get_legend_handles_labels()
    plt.legend(handles = handles, labels = labels)
    fig.savefig('figures/sub/benchmark_' + covariate_name + '.pdf', bbox_inches = 'tight')

