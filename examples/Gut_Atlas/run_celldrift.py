import os
import sys
import random
import pandas as pd 
import numpy as np
import scanpy as sc

def collect_random_cells(adata, key_type = 'cell_type_group', key_perturb = 'Organ', key_time = 'PCW', min_n_cells = 30, normal_n_cells = 100, random_seed = 41):
    df_cells = adata.obs.copy()
    df_stat = pd.DataFrame(df_cells.groupby([key_type, key_time, key_perturb]).size())
    df_stat.columns = ['n_cells']

    cells = []
    for idx_name, n_cells in zip(df_stat.index.values, df_stat['n_cells']):
        name_type, name_time, name_organ = idx_name[0], idx_name[1], idx_name[2]
        smallest_n_across_organs = np.min(df_stat.loc[(name_type, name_time), 'n_cells'])

        # require minimum number of cells in each celltime + time combination
        if smallest_n_across_organs > min_n_cells:
            df_cells_sub = df_cells.loc[(df_cells[key_type] == name_type) & (df_cells[key_perturb] == name_organ) & (df_cells[key_time] == name_time), :].copy()
            random.seed(41)
            cells_sub = random.sample(list(df_cells_sub.index.values), min(n_cells, normal_n_cells)) # get as much cells as possible
            cells += cells_sub 

    adata_sub = adata[cells, :].copy()
    return adata_sub

def run_celldrift(adata, key_type = 'cell_type_group', key_perturb = 'Organ', key_time = 'PCW'):

    sys.path.insert(1, "../test_code/celldrift_code")
    import CellDrift_glm as ct

    adata.obs['size_factor'] = np.sum(adata.X, axis = 1)
    adata.obs['batch'] = 0

    for timepoint in np.unique(adata.obs[key_time]):
        adata_time = adata[adata.obs[key_time] == timepoint, :].copy()
        # run glm model
        adata_time = ct.setup_celldrift(adata_time, cell_type_key = key_type, perturb_key = key_perturb, control_name = 'duojejunum', 
                                    perturb_name = None, size_factor_key = 'size_factor', batch_key = 'batch', min_cells_perGene = 20)
        # compute global coefficients
        suffix_name =  'timepoint_' + timepoint
        adata_time = ct.model_genes(adata_time, n_processes = 16, chunksize = 100, pairwise_contrast_only = True, adjust_batch = False, output_suffix = '_' + suffix_name)
        adata_time = ct.model_selection(adata_time)
        adata_time.write('output_celldrift/' + suffix_name + '.h5ad')

if __name__ == '__main__':
    # load whole file
    adata_all = sc.read('../test_examples/gut_development/data/fetal_RAWCOUNTS_cellxgene.h5ad')
    # get subset
    adata = collect_random_cells(adata_all, key_type = 'cell_type_group', key_perturb = 'Organ', key_time = 'PCW', min_n_cells = 30, random_seed = 41)
    adata.write('fetal_subset_celldrift.h5ad')
    # run celldrift
    run_celldrift(adata)