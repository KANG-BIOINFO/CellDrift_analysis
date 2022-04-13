import pandas as pd
import numpy as np
import scanpy as sc
from anndata import AnnData
import sys
import os

sys.path.insert(1, "../test_code/celldrift_code")
import CellDrift_glm as ct

# list of simulated data
rep_dirs = [i for i in os.listdir("../test_deProb_h5ad_v1/") if os.path.isdir("../test_deProb_h5ad_v1/" + i)]

for rep_dir in rep_dirs:

    h5ad_files = [i for i in os.listdir("../test_deProb_h5ad_v1/" + rep_dir + "/") if i.endswith(".h5ad")]

    # load data and set up celldrift object
    for h5ad_file in h5ad_files:
        DE_size = h5ad_file.split("sim_deProb_")[1].split(".h5ad")[0]
        odata = sc.read("../test_deProb_h5ad_v1/" + rep_dir + "/" + h5ad_file)
        adata = AnnData(X = odata.raw.X.copy(), obs = odata.obs.copy(), var = odata.raw.var.copy())
        adata.obs["Type"] = "Type1"
        adata.obs["Group"] = ["Group" + str(i) for i in adata.obs["Group"]]

        adata = ct.setup_celldrift(adata, cell_type_key = "Type", perturb_key = "Group", control_name = "Group0", 
                                    perturb_name = "Group1", size_factor_key = "nCount", batch_key = "Batch", min_cells_perGene = 5)

        # compute global coefficients
        adata = ct.model_genes(adata, n_processes = 15)
        adata = ct.model_selection(adata)
        adata.write("ct_results/ct_output_DE_" + str(DE_size) + ".h5ad")


    ct_files = [i for i in os.listdir("./ct_results/") if i.endswith(".h5ad")]

    df_pairwise = pd.DataFrame()
    for ct_file in ct_files:
        adata = sc.read("./ct_results/" + ct_file)
        DE_facloc = str(ct_file.split("ct_output_DE_")[1].split(".h5ad")[0])

        df_sub = pd.DataFrame()
        df_sub["contrasts"] = ["Type1_Group1-Type1_Group0"] * adata.shape[1]
        df_sub["mean"] = adata.varm["multcomp_mu"]
        df_sub["lci"] = adata.varm["multcomp_lci"]
        df_sub["uci"] = adata.varm["multcomp_uci"]
        df_sub["SE"] = adata.varm["multcomp_SE"]
        df_sub["p_fdr"] = adata.varm["multcomp_pvals"]
        df_sub["gene"] = adata.var.index.values

        df_sub["DE_FacLoc"] = DE_facloc

        df_pairwise = pd.concat([df_pairwise, df_sub], axis = 0)
    
    df_pairwise.to_csv("ct_results/celldrift_all_Group1_vs_Group0_" + rep_dir + ".txt", sep = "\t")


