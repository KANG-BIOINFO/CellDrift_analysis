import os
import sys
import random
import pandas as pd 
import numpy as np
import scanpy as sc
from tqdm import tqdm

def buildgraph(pairwise_path, suffix):
    output_files = [i for i in os.listdir(pairwise_path) if (i.startswith('glm_predictions_pairwise_comparisons_') and (suffix in i))]

    df_meta = pd.DataFrame(columns = ['time', 'perturbation', 'cell_type'])
    for idx, output_file in enumerate(tqdm(output_files)):
        # get basic information
        df_pairwise = pd.read_csv(pairwise_path + output_file, sep = '\t', header = 0)
        name = output_file.split('glm_predictions_pairwise_comparisons_')[1].split('.txt')[0]
        
        time = name.split('_')[1]

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
    df_combined.to_csv('glm_zscore/pairwise_zscores_combined_' + suffix + '.txt', sep = '\t')
    df_meta.to_csv('glm_zscore/pairwise_contrasts_metadata_' + suffix + '.txt', sep = '\t')


def run_celldrift(adata, suffix_add):

    sys.path.insert(1, "../CellDrift/test_code/celldrift_code")
    import CellDrift_glm as ct
    adata.obs['size_factor'] = np.sum(adata.X, axis = 1)
    adata.obs['batch'] = 0

    for timepoint in np.unique(adata.obs['time']):
        adata_time = adata[adata.obs['time'] == timepoint, :].copy()
        # run glm model
        adata_time = ct.setup_celldrift(adata_time, cell_type_key = 'cell_type', perturb_key = 'perturb', control_name = 'Control', 
                                    perturb_name = None, size_factor_key = 'size_factor', batch_key = 'batch', min_cells_perGene = 0)
        # compute global coefficients
        suffix_name =  'timepoint_' + str(np.round(timepoint, 2)) + '_' + suffix_add
        adata_time = ct.model_genes(adata_time, n_processes = 1, chunksize=60, pairwise_contrast_only = True, adjust_batch = False, output_suffix = '_' + suffix_name)
        adata_time = ct.model_selection(adata_time)
        adata_time.write('output_celldrift/' + suffix_name + '.h5ad')

if __name__ == '__main__':

    h5ad_files = [i for i in os.listdir('./') if i.endswith('.h5ad')]
    for h5ad_file in h5ad_files:
        suffix = h5ad_file.split('.h5ad')[0].split('simulation_')[1]
        # n_times = int(suffix.split('n_times_')[1])

        # load data
        adata = sc.read(h5ad_file)

        # run celldrift
        run_celldrift(adata, suffix_add = suffix)

        # collect zscores
        buildgraph(pairwise_path = './output_celldrift/', suffix = suffix)




