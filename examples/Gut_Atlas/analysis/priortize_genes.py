import pandas as pd 
import numpy as np 

# load data
zscore_file = 'pairwise_zscores_combined_gut_development.txt'
metadata_file = 'pairwise_contrasts_metadata_gut_development.txt'

df_zscores = pd.read_csv(zscore_file, sep = '\t', header = 0, index_col = 0)
df_meta = pd.read_csv(metadata_file, sep = '\t', header = 0, index_col = 0)
df_zscores.values[np.isnan(df_zscores.values)] = 0

# select genes
n = 10
df_genes = pd.DataFrame(columns = ['comparison','time','cell_type','perturbation','genes','sign'])
for comp, time, cell_type, perturbation  in zip(df_meta.index.values, 
                                                df_meta['time'],
                                                df_meta['cell_type'], 
                                                df_meta['perturbation']):

    df_sorted = df_zscores.sort_values(by = [comp], ascending = True).copy()
    down_genes = df_sorted.head(n).index.values
    up_genes = df_sorted.tail(n).index.values 

    df_genes_sub = pd.DataFrame(columns = ['comparison', 'time','perturbation','cell_type','genes','sign'])
    df_genes_sub['comparison'] = [comp] * n * 2
    df_genes_sub['time'] = [time] * n * 2
    df_genes_sub['cell_type'] = [cell_type] * n * 2
    df_genes_sub['perturbation'] = [perturbation] * n * 2
    df_genes_sub['genes'] = list(up_genes) + list(down_genes)
    df_genes_sub['sign'] = ['up'] * n + ['down'] * n 

    df_genes = pd.concat([df_genes, df_genes_sub], axis = 0)

df_genes = df_genes.sort_values(by = ['cell_type', 'perturbation', 'time'], ascending = True)
df_genes.to_csv('Gut_Developement_pairwise_contrast_top_' + str(n) + '.txt', sep = '\t')

