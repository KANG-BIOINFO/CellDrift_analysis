import os
import logging
import pandas as pd
import numpy as np
import scanpy as sc

logging.basicConfig(format='Date-Time : %(asctime)s : Line No. : %(lineno)d - %(message)s', \
                    level = logging.INFO, filename = 'DEG.log')

def format_DEGs(adata):
    keys = ["names","scores","logfoldchanges","pvals","pvals_adj","pts"]
    for i,key in enumerate(keys):
        a = pd.DataFrame(adata.uns["rank_genes_groups"][key]) # transfer to data frame
        b = pd.DataFrame(a.values.T.reshape(1,a.shape[0]*a.shape[1]).T) # reformat the data frame to one column
           
        if i == 0:
            b.columns = [key] # rename the column name
            b["Status"] = sorted(list(a.columns)*a.shape[0]) # add Status annotation
            b.set_index([key],inplace=True)
            b_merged = b
        else:
            if key in ["pts"]:
                pts_all = []
                for cell_group in np.unique(b_merged["Status"]):
                    genes = b_merged.loc[b_merged["Status"] == cell_group,:].index.values
                    pts_all = pts_all + list(a.loc[genes, cell_group])
                b_merged[key] = pts_all
            else:
                b_merged[key] = list(b[0])
        
    return b_merged


if __name__ == '__main__':
    # load data
    if not os.path.isdir('DEGs'):
        os.mkdir('DEGs')

    folders = [
        'simulation_missing_rate_',
        'simulation_n_times',
        'simulation_noise_level',
        'simulation_pattern_gap_',
        'simulation_sequencing_depth_',
        '../simulation_nonLinear/simulation_nonlinear_covariates'
    ]
    covariate_names = ['missing_rate', 'n_times', 'noise', 'pattern_gap', 'sequencing_depth', 'nonlinear_covariate_level']
    methods = ['t-test', 'wilcoxon']

    df_DEG_all = pd.DataFrame()

    for folder, covariate_name in zip(folders, covariate_names):
        logging.info('Run for covariate: ' + covariate_name)
        # load data
        adata_files = [i for i in os.listdir(folder + '/') if i.endswith('.h5ad')]

        for method in methods:
            for adata_file in adata_files:
                logging.info('Run for anndata: ' + adata_file)
                # initialization
                adata = sc.read(folder + '/' + adata_file)
                rep_name = adata_file.split('.h5ad')[0].split('_')[-1]
                covariate_value = adata_file.split('_rep')[0].split(covariate_name)[1]
                # normalize data
                sc.pp.normalize_total(adata, target_sum = 1e6)
                sc.pp.log1p(adata, base = 2)

                # 
                cell_types = adata.obs['cell_type'].unique()
                perturbs = np.setdiff1d(adata.obs['perturb'].unique(), ['Control'])
                timepoints = adata.obs['time'].unique()
                groupby = 'perturb'

                # run though each condition
                for timepoint in timepoints:
                    # print(timepoint)
                    for cell_type in cell_types:
                        # print(cell_type)
                        adata_sub = adata[(adata.obs['time'] == timepoint) & (adata.obs['cell_type'] == cell_type), :].copy()
                        for perturb in perturbs:
                            # print([perturb])
                            # print(adata_sub.shape)
                            sc.tl.rank_genes_groups(adata_sub, groupby = groupby, groups = [perturb], reference = "Control",
                                                    method = method, n_genes = adata.shape[1], pts = True, use_raw = False)
                            df_DEG = format_DEGs(adata_sub)
                            df_DEG['time'] = timepoint
                            df_DEG['cell_type'] = cell_type
                            df_DEG['perturb'] = perturb
                            df_DEG['covariate_name'] = covariate_name
                            df_DEG['covariate_level'] = covariate_value
                            df_DEG['rep'] = rep_name 
                            df_DEG['method'] = method
                            df_DEG_all = pd.concat([df_DEG_all ,df_DEG], axis = 0)

            # write file
    df_DEG_all.to_csv('DEGs/DEG_concat.txt', sep = '\t')
                





