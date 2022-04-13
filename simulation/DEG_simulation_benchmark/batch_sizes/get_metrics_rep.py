import pandas as pd
import numpy as np
import scanpy as sc
import os
from tqdm import tqdm

def TPR(true_up, true_down, pred_up, pred_down):
    # True positive rate or sensivity
    correct_pred = np.union1d((np.intersect1d(true_up, pred_up)) , (np.intersect1d(true_down, pred_down)))
    return len(correct_pred) / len(np.union1d(true_up, true_down))

def FDR(true_up, true_down, true_neg, pred_up, pred_down):
    # False discovery or type I error
    false_pred = np.union1d((np.intersect1d(true_neg, pred_up)) , (np.intersect1d(true_neg, pred_down)))
    correct_pred = np.union1d((np.intersect1d(true_up, pred_up)) , (np.intersect1d(true_down, pred_down)))
    return len(false_pred) / (len(false_pred) + len(correct_pred))

def define_golden(adata, thred):
    adata.var["DEFacGroup_gap"] = [(j-i) for (i,j) in zip(adata.var["DEFacGroup1"], adata.var["DEFacGroup2"])]
    adata.var["abs_DEFacGroup_gap"] = abs(adata.var["DEFacGroup_gap"])

    df_features = adata.var.copy()
    df_features = df_features.sort_values(by = ["DEFacGroup_gap"], ascending = False)
    df_features_up = df_features.loc[df_features["DEFacGroup_gap"] > 0, :]
    df_features_down = df_features.loc[df_features["DEFacGroup_gap"] < 0, :]
    df_features_neg = df_features.loc[df_features["DEFacGroup_gap"] == 0, :]

    up_genes = df_features_up.loc[df_features_up["DEFacGroup_gap"] > np.percentile(df_features_up["DEFacGroup_gap"], 100 - thred), : ].index.values
    down_genes = df_features_down.loc[df_features_down["DEFacGroup_gap"] < np.percentile(df_features_down["DEFacGroup_gap"], thred), : ].index.values
    neg_genes = df_features_neg.loc[(df_features_neg["DEFacGroup1"] == 1) & (df_features_neg["DEFacGroup2"] == 1), : ].index.values

    return up_genes, down_genes, neg_genes

def extract_from_prediction(df_deg, batch_factor, category, p_fdr_thred = 0.05):
    df_sub = df_deg.loc[df_deg["Batch_FacLoc"] == batch_factor, : ].copy()
    
    if category == "pairwise":
        up_genes = df_sub.loc[(df_sub["p_fdr"] < p_fdr_thred) & (df_sub["mean"] > 0), "gene"]
        down_genes = df_sub.loc[(df_sub["p_fdr"] < p_fdr_thred) & (df_sub["mean"] < 0), "gene"]
        neg_genes = df_sub.loc[(df_sub["p_fdr"] > p_fdr_thred), "gene"]

    elif category == "others":
        up_genes = df_sub.loc[(df_sub["pvals_adj"] < p_fdr_thred) & (df_sub["scores"] > 0), :].index.values
        down_genes = df_sub.loc[(df_sub["pvals_adj"] < p_fdr_thred) & (df_sub["scores"] < 0), :].index.values
        neg_genes = df_sub.loc[(df_sub["pvals_adj"] > p_fdr_thred), :].index.values
    
    elif category == "MAST":
        up_genes = df_sub.loc[(df_sub["Pr(>Chisq)"] < p_fdr_thred) & (df_sub["coef"] > 0), "primerid"]
        down_genes = df_sub.loc[(df_sub["Pr(>Chisq)"] < p_fdr_thred) & (df_sub["coef"] < 0), "primerid"]
        neg_genes = df_sub.loc[(df_sub["Pr(>Chisq)"] > p_fdr_thred), "primerid"]

    return up_genes, down_genes, neg_genes


if __name__ == "__main__":

    # load data
    threshold = 75
    raw_path = "../test/test_benchmark/simulated_ifn_data/test_batch_h5ad_v2/"
    
    df_output = pd.DataFrame(columns = ["TPR", "FDR", "Replicates", "BatchFacLoc", "Method"])

    row_id = 0
    for rep_id in tqdm(np.arange(1,11)):
        rep_name = "Rep" + str(rep_id)
        df_pairwise = pd.read_csv("ct_results/celldrift_all_Group1_vs_Group0_" + rep_name + ".txt", sep = "\t", header = 0, index_col = 0)
        df_t = pd.read_csv("others_results/" + rep_name + "_t-test_all_Group1_vs_Group0.txt", sep = "\t", header = 0, index_col = 0)
        df_wilcox = pd.read_csv("others_results/" + rep_name + "_wilcoxon_all_Group1_vs_Group0.txt", sep = "\t", header = 0, index_col = 0)
        df_t_over = pd.read_csv("others_results/" + rep_name + "_t-test_overestim_var_all_Group1_vs_Group0.txt", sep = "\t", header = 0, index_col = 0)
        df_mast = pd.read_csv("mast_results/" + rep_name + "_mast_Group2_vs_Group1.txt", sep = "\t", header = 0, index_col = 0)
        
        df_others = [df_t, df_wilcox, df_t_over]
        test_names = ["t-test", "wilcoxon", "t-test-overdispersion"]

        ct_files = [i for i in os.listdir(raw_path + rep_name + "/") if i.endswith(".h5ad")]
        
        # calculate metrics
        for ct_file in ct_files:
            adata = sc.read(raw_path + rep_name + "/" + ct_file)
            batch_facloc = float(ct_file.split("sim_batchLoc_")[1].split(".h5ad")[0])
            
            genes_up, genes_down, genes_neg = define_golden(adata, thred = threshold)

            # ct pairwise
            pred_up, pred_down, pred_neg = extract_from_prediction(df_pairwise, batch_facloc, "pairwise")
            TPR_score = TPR(genes_up, genes_down, pred_up, pred_down)
            FDR_score = FDR(genes_up, genes_down, genes_neg, pred_up, pred_down)
            
            df_output.loc[row_id, :] = [TPR_score, FDR_score, rep_name, batch_facloc, "CellDrift"]
            row_id += 1

            # others
            for name, df_deg in zip(test_names, df_others):
                pred_up, pred_down, pred_neg = extract_from_prediction(df_deg, batch_facloc, "others")
                TPR_score = TPR(genes_up, genes_down, pred_up, pred_down)
                FDR_score = FDR(genes_up, genes_down, genes_neg, pred_up, pred_down)

                df_output.loc[row_id, :] = [TPR_score, FDR_score, rep_name, batch_facloc, name]
                row_id += 1
            
            # mast
            pred_up, pred_down, pred_neg = extract_from_prediction(df_mast, batch_facloc, "MAST")
            TPR_score = TPR(genes_up, genes_down, pred_up, pred_down)
            FDR_score = FDR(genes_up, genes_down, genes_neg, pred_up, pred_down)
            df_output.loc[row_id, :] = [TPR_score, FDR_score, rep_name, batch_facloc, "MAST"]
            row_id += 1

    df_output.to_csv("BenchBatches_metrics_threshold_" + str(threshold) + ".txt", sep = "\t")

        

