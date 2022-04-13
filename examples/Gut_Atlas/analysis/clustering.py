import os
import pickle
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import adjusted_rand_score

import skfda
from skfda import datasets
from skfda.exploratory.visualization.clustering import (
    ClusterMembershipLinesPlot,
    ClusterMembershipPlot,
    ClusterPlot,
)
from skfda.ml.clustering import FuzzyCMeans, KMeans

colormap = plt.cm.get_cmap('tab20b')

# create a fd object
def draw_eachCluster(fd, labels, cluster_id = 0, n_clusters = 20, prefix_name = 'figures_cluster_ind_'):
    labels = pd.Series(labels, dtype = 'category')
    indices_samples = np.arange(len(labels))[labels == cluster_id]
    
    fd_cluster = fd[indices_samples]
    labels_cluster = labels[indices_samples].cat.remove_unused_categories()

    cluster_colors = colormap(np.arange(n_clusters))

    fig, ax = plt.subplots()
    fd_cluster.plot(group = np.array(labels_cluster.values), 
                    group_names = labels_cluster.cat.categories,
                    group_colors = cluster_colors,
                    axes = ax)
    fig.savefig('figures/fuzzy_clusters/' + prefix_name + str(cluster_id), bbox_inches = 'tight')
    plt.close()

def perform_clustering(df_zscores, df_metadata, genes, cell_type = 'cMono', 
                        perturbation = 'COVID_SEV', n_clusters = 20, seed = 42,
                        fig_name = 'figures_cluster'):
    print(cell_type)
    print(perturbation)

    # preparation
    df_metadata['time_v2'] = [float(i.split('F')[1]) for i in df_metadata['time']]
    df_metadata = df_metadata.loc[(df_metadata['cell_type'] == cell_type) & (df_metadata['perturbation'] == perturbation), :].copy()
    df_metadata = df_metadata.sort_values(['time_v2'], ascending = True)
    timepoints = list(df_metadata['time_v2'])

    comps = list(df_metadata.index.values)
    df_zscores = df_zscores[comps]
    df_zscores = df_zscores.loc[genes, :]

    data_matrix = []
    for gene in df_zscores.index.values:
        gene_zvalues = df_zscores.loc[gene, :]
        data_matrix.append(gene_zvalues)

    # create representation of functional data
    fd = skfda.FDataGrid(
        data_matrix = data_matrix,
        grid_points = timepoints
    )
    
    # run clustering methods
    kmeans = KMeans(n_clusters=n_clusters, random_state=seed)
    kmeans.fit(fd)
    pred_labels = kmeans.predict(fd)
    
    fuzzy_kmeans = FuzzyCMeans(n_clusters=n_clusters, random_state=seed)
    fuzzy_kmeans.fit(fd)
    pred_labels_fuzzy = fuzzy_kmeans.predict(fd)

    # visualizations
    fig, ax = plt.subplots()
    ClusterPlot(fuzzy_kmeans, fd, axes = ax).plot()
    fig.savefig('figures/fuzzy_clusters/' + fig_name + '.pdf', bbox_inches = 'tight')

    # for cluster_id in range(n_clusters):
        # draw_eachCluster(fd, pred_labels_fuzzy, cluster_id = cluster_id, prefix_name = 'figures_cluster_ind_')

    df_clusters = pd.DataFrame()
    df_clusters['genes'] = genes
    df_clusters['clusters_kmeans'] = pred_labels
    df_clusters['clusters_fuzzy'] = pred_labels_fuzzy

    return fd, df_clusters

if __name__ == '__main__':
    # initialization
    df_zscores = pd.read_csv('pairwise_zscores_combined_gut_development.txt', sep = '\t', header = 0, index_col = 0)
    df_metadata = pd.read_csv('pairwise_contrasts_metadata_gut_development.txt', sep = '\t', header = 0, index_col = 0)
    df_features = pd.read_csv('significant_genes.txt', sep = '\t', header = 0, index_col = 0)
    df_zscores.values[np.isnan(df_zscores.values)] = 0

    # run clustering
    cell_types = ['epithelium', 'mesenchymal']
    perturbs = ['colon', 'ileum']
    n_clusters = 20

    for cell_type in cell_types:
        for perturb in perturbs:
            group = [i for i in df_features['group'].unique() if (cell_type in i) and (perturb in i)][0]
            features = df_features.loc[df_features['group'] == group, 'genes']
            fd, df_clusters = perform_clustering(df_zscores, df_metadata, features,
                                                cell_type = cell_type, perturbation = perturb, 
                                                n_clusters = n_clusters, seed = 42)
            df_clusters.to_csv('fuzzy_clusters/clusters_' + cell_type + '_' + perturb + '.txt', sep = '\t')


