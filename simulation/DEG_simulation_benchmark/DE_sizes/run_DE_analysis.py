# test DEG approaches (wilcoxon, t-test, t-test-overdisp) on simulated data
import scanpy as sc
import numpy as np
import pandas as pd
from anndata import AnnData
import os

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

### run
source_path = "../test/test_benchmark/simulated_ifn_data/test_deProb_h5ad_v1/"
rep_dirs = [i for i in os.listdir(source_path) if os.path.isdir(source_path + i)]

for rep_dir in rep_dirs:

    h5ad_files = [i for i in os.listdir(source_path + rep_dir + "/") if i.endswith(".h5ad")]
    methods = ["t-test", "wilcoxon", "t-test_overestim_var"] # logreg

    for method in methods:
        df_DEG_all = pd.DataFrame()
        
        for h5ad_file in h5ad_files:
            DE_facLoc = h5ad_file.split("sim_deProb_")[1].split(".h5ad")[0]

            odata = sc.read(source_path + rep_dir + "/" + h5ad_file)
            adata = AnnData(X = odata.raw.X.copy(), obs = odata.obs.copy(), var = odata.raw.var.copy())

            adata.obs["Group"] = ["Group" + str(i) for i in adata.obs["Group"]]

            sc.pp.normalize_total(adata, target_sum = 1e6)
            sc.pp.log1p(adata, base = 2)

            sc.tl.rank_genes_groups(adata, "Group", groups = ["Group1"], reference = "Group0", method = method, n_genes = adata.shape[1], pts = True, use_raw = False)
            df_DEG = format_DEGs(adata)
            df_DEG["DE_FacLoc"] = DE_facLoc
            df_DEG_all = pd.concat([df_DEG_all, df_DEG], axis = 0)
            
        df_DEG_all.to_csv("others_results/" + rep_dir + "_" + method + "_all_Group1_vs_Group0.txt", sep = "\t")