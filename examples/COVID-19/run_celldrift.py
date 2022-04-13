import os
import sys
import random
import pandas as pd 
import numpy as np
import scanpy as sc

def run_celldrift(adata_all_path:

	sys.path.insert(1, "../test_code/celldrift_code")
	import CellDrift_glm as ct

	# load data
	adata_all = sc.read(adata_all_path)
	# sanity check
	print(np.sum(adata_all.X[0,:].toarray()))

	# instead of using variable genes, we use all genes
	bdata = sc.read('../../glm/covid_atlas_all_genes.h5ad')
	genes = [i for i in bdata.var.index.values if (not i.startswith('AB_'))]
	genes = np.intersect1d(genes, adata_all.var.index.values)
	print(len(genes))

	adata_all = adata_all[:, genes].copy()
	sc.pp.filter_genes(adata_all, min_cells = 0)
	print(adata_all.shape)

	adata_ctrl = adata_all[adata_all.obs['Source'] == 'HV', : ].copy()
	adata_pert = adata_all[adata_all.obs['Source'] != 'HV', : ].copy()

	for timepoint in np.unique(adata_pert.obs['TimeSinceOnset']):
		adata_time = adata_pert[adata_pert.obs['TimeSinceOnset'] == timepoint, :].copy()
		adata_time.obs['Source'] = adata_time.obs['Source'].cat.remove_unused_categories()
		for pert in np.unique(adata_time.obs['Source']):
			print(timepoint, pert)
			adata = adata_time[adata_time.obs['Source'] == pert, :].copy().concatenate(adata_ctrl.copy(), join = 'outer', index_unique = None)
			shared_types = [i for i in np.unique(adata.obs['Annotation_major_subset']) if (pd.DataFrame(adata.obs.groupby(['Annotation_major_subset', 'Source']).size()).loc[i,:].shape[0] == 2)]
			adata = adata[adata.obs['Annotation_major_subset'].isin(shared_types), :].copy()

			# run glm model
			adata = ct.setup_celldrift(adata, cell_type_key = "Annotation_major_subset", perturb_key = 'Source', control_name = "HV", 
										perturb_name = None, size_factor_key = "QC_total_UMI", batch_key = "scRNASeq_sample_ID", min_cells_perGene = 10)
			# compute global coefficients
			suffix_name =  '_d' + str(int(float(timepoint))) + '_' + pert
			adata = ct.model_genes(adata, n_processes = 16, pairwise_contrast_only = True, adjust_batch = False, output_suffix = suffix_name)
			adata = ct.model_selection(adata)
			adata.write('covid_atlas_all_HVGs' + suffix_name + '.h5ad')

if __name__ == '__main__':
	run_celldrift('COMBAT-CITESeq-DATA_selection_v2_time_raw_rep3.h5ad')