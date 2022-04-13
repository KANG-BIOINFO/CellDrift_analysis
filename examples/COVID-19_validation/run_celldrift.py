import os
import sys
import random
import pandas as pd 
import numpy as np
import scanpy as sc

def run_celldrift(adata_all_path, rep_name):

	# create new folder
	sys.path.insert(1, "/Volumes/aronow/Kang/CellDrift/test_code/celldrift_code")
	import CellDrift_glm as ct

	# load data
	adata_all = sc.read(adata_all_path)

	# QC
	sc.pp.filter_genes(adata_all, min_cells = 100)
	sc.pp.filter_cells(adata_all, min_counts = 0)
	adata_all.obs['batch'] = 0
	print(adata_all.shape)

	# get control
	adata_ctrl = adata_all[adata_all.obs['condition'] == 'control', :].copy()

	perturbs = np.setdiff1d(adata_all.obs['condition'].unique(), ['control'])
	for perturb in perturbs:
		adata_pert = adata_all[adata_all.obs['condition'] == perturb, :].copy()
		timepoints = adata_pert.obs['time'].unique()
		for timepoint in timepoints:
			# get concatenated input for celldrift
			adata_time = adata_pert[adata_pert.obs['time'] == timepoint, :].copy()
			adata_input = adata_ctrl.concatenate(adata_time, join = 'outer', index_unique = None)
			# adata_input.obs['condition'] = adata_input.obs['condition'].cat.remove_unused_categories()

			adata_input = ct.setup_celldrift(adata_input, cell_type_key = "cell_type", perturb_key = 'condition', control_name = "control", 
							perturb_name = None, size_factor_key = "n_counts", batch_key = "batch", min_cells_perGene = 10)
			# compute global coefficients
			suffix_name =  '_d' + str(timepoint) + '_' + perturb + '_' + rep_name
			adata_input = ct.model_genes(adata_input, n_processes = 16, chunksize = 100, pairwise_contrast_only = True, adjust_batch = False, output_suffix = suffix_name)
			adata_input = ct.model_selection(adata_input)
			adata_input.write('output_celldrift/covid_Ren' + suffix_name + '.h5ad')

if __name__ == '__main__':
	rep_names = ['rep0', 'rep1', 'rep2']		
	rep_path_list = ['../data/Ren_large_covidatlas_raw_selection_rep0.h5ad',
					'../data/Ren_large_covidatlas_raw_selection_rep1.h5ad',
					'../data/Ren_large_covidatlas_raw_selection_rep2.h5ad']
	for rep_name, rep_path in zip(rep_names, rep_path_list):
		run_celldrift(rep_path, rep_name)