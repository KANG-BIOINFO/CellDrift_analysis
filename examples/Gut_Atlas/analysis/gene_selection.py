import os
import pickle
import scanpy as sc
import numpy as np 
import pandas as pd 
import itertools

p_thred = 0.05
pairwise_files = [i for i in os.listdir('../output_celldrift/') if i.startswith('glm_predictions_pairwise')]

genes_significant = {}
for pairwise_file in pairwise_files:
    df_pairwise = pd.read_csv('../output_celldrift/' + pairwise_file, sep = '\t', header = 0, index_col = 0)
    df_pairwise['cell_type_v2'] = [(i.split('\',')[0].split('(\'')[1]) for i in df_pairwise['cell_type']]
    df_pairwise['perturbation_v2'] = [(i.split('\',')[0].split('(\'')[1]) for i in df_pairwise['perturbation']]
    
    cell_types = list(df_pairwise['cell_type_v2'].unique())
    perturbations = list(df_pairwise['perturbation_v2'].unique())
    groups = list(itertools.product(cell_types, perturbations))
    for group in groups:
        df_sub = df_pairwise.loc[(df_pairwise['cell_type_v2']==group[0]) & (df_pairwise['perturbation_v2'] == group[1]), :].copy()
        sign_genes = df_sub.loc[df_sub['p_fdr'] < p_thred, 'gene']
        if group not in genes_significant.keys():
            genes_significant[group] = sign_genes
        else:
            genes_significant[group] = np.union1d(sign_genes, genes_significant[group])

genes = []; gene_groups = []
for group in genes_significant.keys():
    print(group)
    print(len(genes_significant[group]))
    genes += list(genes_significant[group])
    gene_groups += [group] * len(genes_significant[group])

df_genes = pd.DataFrame(columns = ['genes', 'group'])
df_genes['genes'] = genes 
df_genes['group'] = gene_groups
df_genes.to_csv('significant_genes.txt', sep = '\t')