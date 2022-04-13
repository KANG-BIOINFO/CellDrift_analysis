import os
import sys
import random
import pandas as pd 
import numpy as np
import scanpy as sc
from tqdm import tqdm

def buildgraph(pairwise_path, suffix):
    output_files = [i for i in os.listdir(pairwise_path) if (i.startswith('glm_predictions_pairwise_comparisons_'))]

    df_meta = pd.DataFrame(columns = ['time', 'perturbation', 'cell_type'])
    for idx, output_file in enumerate(tqdm(output_files)):
        # get basic information
        df_pairwise = pd.read_csv(pairwise_path + output_file, sep = '\t', header = 0)
        time = output_file.split('glm_predictions_pairwise_comparisons_timepoint_')[1].split('.txt')[0]
        
        n_types = len(np.unique(df_pairwise['cell_type']))
        n_perts = len(np.unique(df_pairwise['perturbation']))
        n_genes = len(np.unique(df_pairwise['gene']))

        # make graph per (time + pert)
        df_pairwise['zscore'] = df_pairwise['z']
        coeff_vals = df_pairwise['zscore'].values.reshape([n_genes, n_types * n_perts])

        genes = list(df_pairwise.sort_values(by = ['cell_type', 'perturbation', 'gene'])['gene'][ : n_genes])
        contrasts = [(time + '_' + i) for i in list(df_pairwise['contrast'][ : n_types * n_perts])]
        avail_types = [(i.split('\',')[0].split('(\'')[1]) for i in list(df_pairwise['cell_type'][:(n_types * n_perts)])]
        avail_perts = [(i.split('\',')[0].split('(\'')[1]) for i in list(df_pairwise['perturbation'][:(n_types * n_perts)])]
        
        df_graph = pd.DataFrame(data = coeff_vals, index = genes, columns = contrasts)

        if idx == 0:
            df_combined = df_graph.copy()
        else:
            df_combined = pd.concat([df_combined, df_graph], axis = 1, join = 'outer')
        
        print(contrasts)
        for contrast, cell_type, pert in zip(contrasts, avail_types, avail_perts):
            df_meta.loc[contrast, : ] =  [time, pert, cell_type]
    
    if not os.path.isdir('glm_zscore'):
        os.mkdir('glm_zscore')

    df_combined.to_csv('pairwise_zscores_combined_' + suffix + '.txt', sep = '\t')
    df_meta.to_csv('pairwise_contrasts_metadata_' + suffix + '.txt', sep = '\t')

if __name__ == '__main__':
    pairwise_path = '../output_celldrift/'
    suffix = 'gut_development'
    buildgraph(pairwise_path, suffix)