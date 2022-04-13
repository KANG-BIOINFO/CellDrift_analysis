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

def calculate_scores(df_deg, cell_type = 'Type_0', perturbation = 'Perturb_0',
                    n_clusters = 3, n_genes_pattern = 20, seed = 42):
    
    # preparation of input
    df_deg = df_deg.loc[(df_deg['cell_type'] == cell_type) & (df_deg['perturb'] == perturbation), : ].copy()
    df_deg = df_deg.sort_values(['genes_order', 'time'], ascending = True)

    timepoints = np.unique(df_deg['time'])
    timepoints.sort()

    # initialize a data matrix
    data_matrix = []
    n_genes = n_clusters * n_genes_pattern

    for i in range(n_genes):
        gene_name = 'Gene_' + str(i)
        df_scores = df_deg.loc[df_deg['names'] == gene_name, :]
        data_matrix.append(df_scores['scores'])                                         

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
    df_DEG = pd.read_csv('../DEGs/DEG_concat.txt', sep = '\t', header = 0)
    df_DEG['covariate_level'] = [i.split('_')[1] for i in df_DEG['covariate_level']]
    df_DEG['genes_order'] = [int(i.split('_')[1]) for i in df_DEG['names']]
    
    methods = ['t-test', 'wilcoxon']
    cell_types = ['Type_0', 'Type_1']
    perturbs = ['Perturb_0', 'Perturb_1']    

    covariate_names = ['missing_rate', 'n_times', 'noise', 'pattern_gap', 'sequencing_depth', 'nonlinear_covariate_level']
    columns = ['covariate_name', 'covariate_level', 'rep', 'cell_type', 'perturb', 'ARI', 'correlation', 'method']
    df_metrics = pd.DataFrame(columns = columns)
    n_reps = 3

    idx = 0
    for method in methods:
        print(method)
        df_DEG_sub1 = df_DEG.loc[df_DEG['method'] == method, :].copy()

        for covariate_name in covariate_names:
            print(covariate_name)
            df_DEG_sub2 = df_DEG_sub1.loc[df_DEG_sub1['covariate_name'] == covariate_name, :].copy()
            covariate_values = np.unique(df_DEG_sub2['covariate_level'])

            for covariate_value in covariate_values:
                print(covariate_value)
                df_DEG_sub3 = df_DEG_sub2.loc[df_DEG_sub2['covariate_level'] == covariate_value, : ].copy()

                for cell_type in cell_types:
                    for perturb in perturbs:

                        for rep in range(n_reps):
                            if (covariate_name == 'sequencing_depth') and (rep in [1,2]):
                                continue
                            rep_name = 'rep' + str(rep)
                            df_DEG_sub4 = df_DEG_sub3.loc[df_DEG_sub3['rep'] == rep_name, :].copy()
                            df_DEG_sub4 = df_DEG_sub4.sort_values(by = ['genes_order', 'time'], ascending = True)
                            df_DEG_sub4.to_csv('test.txt', sep = '\t')
                            ARI_score_kmeans, ARI_score_fuzzy = calculate_scores(df_DEG_sub4, cell_type = cell_type, perturbation = perturb)

                            # append results
                            result_kmeans = [covariate_name, covariate_value, rep, cell_type, perturb, ARI_score_kmeans, 0, 'kmeans_DEGs_' + method]
                            df_metrics.loc[idx, :] = result_kmeans
                            idx += 1

                            result_fuzzy = [covariate_name, covariate_value, rep, cell_type, perturb, ARI_score_fuzzy, 0, 'kmeans_fuzzy_DEGs_' + method]
                            df_metrics.loc[idx, :] = result_fuzzy 
                            idx += 1
    
    df_metrics.to_csv('results_clusters_kmeans_fuzzy_DEGs.txt', sep = '\t')


    


