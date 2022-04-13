import os
from tqdm import tqdm
import numpy as np
import pandas as pd 
import scanpy as sc

output_files = [i for i in os.listdir('../output_celldrift') if i.startswith('glm_predictions_pairwise_comparisons_')]

df_meta = pd.DataFrame(columns = ['time', 'perturbation', 'cell_type', 'rep'])
for idx, output_file in enumerate(tqdm(output_files)):
    # get basic information
    df_pairwise = pd.read_csv('../output_celldrift/' + output_file, sep = '\t', header = 0)
    name = output_file.split('glm_predictions_pairwise_comparisons_')[1].split('.txt')[0]
    
    time = name.split('_')[0]
    time_v2 = int(time.split('d')[1])
    rep_name = name.split('_')[-1]
    pert = 'mild' if 'mild' in name else 'severe'

    n_types = len(np.unique(df_pairwise['cell_type']))
    n_genes = len(np.unique(df_pairwise['gene']))

    # make graph per (time + pert)
    df_pairwise['zscore'] = df_pairwise['mean'] / df_pairwise['SE']
    coeff_vals = df_pairwise['zscore'].values.reshape([n_genes, n_types])

    genes = list(df_pairwise.sort_values(by = ['cell_type', 'gene'])['gene'][ : n_genes])
    contrasts = [(time + '_' + i + '_' + rep_name) for i in list(df_pairwise['contrast'][ : n_types])]
    avail_types = [(i.split('\',')[0].split('(\'')[1]) for i in list(df_pairwise['cell_type'][:n_types])]
    
    df_graph = pd.DataFrame(data = coeff_vals, index = genes, columns = contrasts)

    if idx == 0:
        df_combined = df_graph.copy()
    else:
        df_combined = pd.concat([df_combined, df_graph], axis = 1, join = 'outer')
    
    for contrast, cell_type in zip(contrasts, avail_types):
        df_meta.loc[contrast, : ] =  [time, pert, cell_type, rep_name]

print(np.sum(np.sum(np.isnan(df_combined.values))))
df_combined.to_csv('pairwise_zscores_combined.txt', sep = '\t')
df_meta.to_csv('pairwise_contrasts_metadata.txt', sep = '\t')




