import os
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

def calculate_scores(df_zscore, df_meta, cell_type = 'Type_0', perturbation = 'Perturb_0',
                    n_clusters = 3, n_genes_pattern = 20, seed = 42):
    
    # preparation of input
    df_meta = df_meta.loc[(df_meta['cell_type'] == cell_type) & (df_meta['perturbation'] == perturbation), :].copy()
    df_meta = df_meta.sort_values(by = ['time'], ascending = True)

    timepoints = df_meta['time']
    df_zscore = df_zscore[list(df_meta.index.values)]
    df_zscore.values[np.isnan(df_zscore.values)] = 0                        

    # initialize a data matrix
    data_matrix = []
    for i in range(df_zscore.shape[0]):
        gene_name = 'Gene_' + str(i)
        gene_zvalues = df_zscore.loc[gene_name, :]
        data_matrix.append(gene_zvalues)
                        

    # create representation of functional data
    fd = skfda.FDataGrid(
        data_matrix = data_matrix,
        grid_points = timepoints
    )

    # create true labels
    true_labels = np.repeat(np.arange(n_clusters), n_genes_pattern)

    # run kmeans
    kmeans = KMeans(n_clusters=n_clusters, random_state=seed)
    kmeans.fit(fd)
    pred_labels = kmeans.predict(fd)
    ARI_score_kmeans = adjusted_rand_score(true_labels, pred_labels)

    # run fuzzy kmeans
    fuzzy_kmeans = FuzzyCMeans(n_clusters=n_clusters, random_state=seed)
    fuzzy_kmeans.fit(fd)
    pred_labels_fuzzy = fuzzy_kmeans.predict(fd)
    ARI_score_fuzzy = adjusted_rand_score(true_labels, pred_labels_fuzzy)
    return (ARI_score_kmeans, ARI_score_fuzzy)

if __name__ == '__main__':
    # load data
    folders = [
        'simulation_missing_rate_',
        'simulation_n_times',
        'simulation_noise_level',
        'simulation_pattern_gap_',
        'simulation_sequencing_depth_',
        '../simulation_nonLinear/simulation_nonlinear_covariates'
    ]
    covariate_names = ['missing_rate', 'n_times', 'noise', 'pattern_gap', 'sequencing_depth', 'nonlinear_covariate_level']
    
    cell_types = ['Type_0', 'Type_1']
    perturbs = ['Perturb_0', 'Perturb_1']
    columns = ['covariate_name', 'covariate_level', 'rep', 'cell_type', 'perturb', 'ARI', 'correlation', 'method']
    df_metrics = pd.DataFrame(columns = columns)

    idx = 0
    for folder, covariate_name in zip(folders, covariate_names):
        print(covariate_name)

        zscore_files = [i for i in os.listdir('../' + folder + '/glm_zscore/') if i.startswith('pairwise_zscores_combined')]
        meta_files = [i for i in os.listdir('../' + folder + '/glm_zscore/') if i.startswith('pairwise_contrasts_metadata')]

        zscore_files.sort()
        meta_files.sort()
 
        for zscore_file, meta_file in zip(zscore_files, meta_files):
            df_zscore = pd.read_csv('../' + folder + '/glm_zscore/' + zscore_file, sep = '\t', header = 0, index_col = 0)
            df_meta = pd.read_csv('../' + folder + '/glm_zscore/' + meta_file, sep = '\t', header = 0, index_col = 0)

            covariate_value = zscore_file.split('_rep')[0].split(covariate_name + '_')[1]
            rep_id = int(zscore_file.split('_rep')[1].split('.txt')[0])

            for cell_type in cell_types:
                for perturb in perturbs:
                    ARI_kmeans, ARI_fuzzy = calculate_scores(df_zscore, df_meta, cell_type = cell_type, perturbation = perturb,
                                                                n_clusters = 3, n_genes_pattern = 20, seed = 42)
                    
                    result_kmeans = [covariate_name, covariate_value, rep_id, cell_type, perturb, ARI_kmeans, 0, 'kmeans']
                    df_metrics.loc[idx, :] = result_kmeans
                    idx += 1

                    result_fuzzy = [covariate_name, covariate_value, rep_id, cell_type, perturb, ARI_fuzzy, 0, 'kmeans_fuzzy']
                    df_metrics.loc[idx, :] = result_fuzzy 
                    idx += 1
    
    df_metrics.to_csv('results_clusters_kmeans_fuzzy.txt', sep = '\t')


    


