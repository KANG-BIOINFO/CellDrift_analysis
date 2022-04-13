library(stringr)
library(fdapace)
library(EMCluster)
library(mclust)
library(dplyr)
library(ggplot2)
library(cowplot)

setwd('../test_examples/simulation_publication/simulation_v3/add_replicates/')

# preprocess zscore files and metadata
preprocess_input = function(df_combined_path, df_meta_path)
{
  df_combined = read.table(df_combined_path, sep = '\t', header = 1, row.names = 1)
  df_meta = read.table(df_meta_path, sep = '\t', header = 1, row.names = 1)
  
  rownames(df_meta) = paste0('X', str_replace(rownames(df_meta), '-', '.'))
  df_meta$id = paste0(df_meta$cell_type, '_', df_meta$perturbation)
  head(df_meta)
  
  genes = c()
  n_genes = dim(df_combined)[1]
  for (i in 0:(n_genes - 1)){genes = append(genes, paste0('Gene_', i))}
  df_combined = df_combined[genes,]
  
  df_list = list()
  df_list[['df_combined']] = df_combined
  df_list[['df_meta']] = df_meta
  return(df_list)
}

# prepare fpca input data
prepare_fpca = function(df_zscore, df_meta, genes, cell_type = 'Type_0', perturb = 'Perturb_0'){
  # prepare input data 
  df_submeta = df_meta[(df_meta$cell_type == cell_type) & (df_meta$perturbation == perturb), ]
  df_submeta = df_submeta[order(df_submeta$time),]
  
  comparisons = rownames(df_submeta)
  df_zscore = df_zscore[genes, ]
  df_zscore = df_zscore[, comparisons]
  
  yVec = c(); tVec = c(); IDs = c()
  for (i in 1:ncol(df_zscore)){
    col = colnames(df_zscore)[i]
    yVec = append(yVec, df_zscore[,i])
    tVec = append(tVec, rep(df_submeta[col, 'time'], nrow(df_zscore)))
    IDs = append(IDs, rownames(df_zscore))
  }
  input = MakeFPCAInputs(IDs = IDs,
                         tVec = tVec,
                         yVec = yVec)
  return (input)
}

# calculate ARI score for cluster prediction
calculate_cluster_accuracy = function(cout, n_patterns)
{
  n_genes = length(cout$cluster)
  predicted_labels = cout$cluster
  true_labels = c()
  n_genes_pattern = as.integer(n_genes / n_patterns)
  for (i in 1:n_patterns)
  {
    true_labels = append(true_labels, rep(i, n_genes_pattern))
  }
  ari = adjustedRandIndex(true_labels, predicted_labels)
  return(ari)
}

# calculate MSE for imputation errors
calculate_correlation = function(timepoints, fpca_obj, true_vals){
  # initialization
  xiEst = fpca_obj$xiEst
  phi = fpca_obj$phi
  workgrids = fpca_obj$workGrid
  n_grids = length(workgrids)
  mu = fpca_obj$mu
  
  # retrieve values
  for (i in 1:length(timepoints))
  {
    timepoint = timepoints[i]
    # go for the PC score of cloest work grid
    if (timepoint <= workgrids[1])
    {
      pred_val = mu[1] + xiEst %*% phi[1,]
    } 
    else if (timepoint == 1)
    {
      pred_val = mu[n_grids] + xiEst %*% phi[n_grids,]
    } 
    # go for weighted score between two work grid
    else {
      query_index = sum((workgrids - timepoint) <= 0)
      if (query_index == length(workgrids)){break}
      g1 = workgrids[query_index]
      g2 = workgrids[query_index + 1]
      gap = g2 - g1
      phi_val = phi[query_index, ] * ((g2 - timepoint) / gap ) + phi[query_index + 1, ] * ((timepoint - g1) / gap)
      mu_val = mu[query_index] * ((g2 - timepoint) / gap ) + mu[query_index + 1] * ((timepoint - g1) / gap)
      pred_val = mu_val + xiEst %*% phi_val
    }
    # combine results
    if (i == 1)
    {
      pred_vals = pred_val
    }
    else{
      pred_vals = cbind(pred_vals, pred_val)
    }
  }
  pred_vals_sub = pred_vals[append(1:20, 41:60), ]
  true_vals_sub = true_vals[append(1:20, 41:60), ]
  true_vals_sub = true_vals_sub[,  1:dim(pred_vals_sub)[2]]
  
  avg_corr = mean(diag(cor(t(pred_vals_sub), t(true_vals_sub), method = 'pearson')))
  return(avg_corr)
}

retrieve_true_labels = function(type_pattern, n_timepoints, n_genes_pattern, mean_pert, mean_type_pert, pattern_coeff_gap )
{
  timepoints = seq(0,1,1/(n_timepoints-1))
  true_vals = t(replicate(n_genes_pattern, mean_pert + pattern_coeff_gap * mean_type_pert - pattern_coeff_gap * mean_type_pert *timepoints))
  true_vals = rbind(true_vals, t(replicate(n_genes_pattern, rep(mean_pert + mean_type_pert, n_timepoints))))
  true_vals = rbind(true_vals, t(replicate(n_genes_pattern, mean_pert + pattern_coeff_gap * mean_type_pert + pattern_coeff_gap * mean_type_pert *timepoints)))
  return (true_vals)
}

run_pipeline_over_zscore = function(covariate_key = 'missing_rate', covariate_values = seq(0.0, 0.9, 0.1),
                                    folder = '../glm_zscore/', pattern_coeff_gap = 1.0, mean_pert = 0.3, 
                                    mean_type_pert = 2.0, n_patterns = 3, n_genes_pattern = 20, n_timepoints = 20, n_rep = 3,
                                    cmethod = 'EMCluster', seed = 42)
{
  # initialization
  df_metrics = data.frame(matrix(ncol = 8, nrow = 0))
  colnames(df_metrics) = c('covariate_name', 'covariate_level', 'rep', 'cell_type', 'perturb', 'ARI', 'correlation', 'method')
  
  # get true time coefficients based on simulation parameters
  for (covariate_val in covariate_values)
  {
    print(covariate_val)
    for (cell_type in c('Type_0','Type_1'))
    {
      for (perturb in c('Perturb_0', 'Perturb_1'))
      {
        for (rep_id in 0:(n_rep-1))
        {
          set.seed(seed)
          print(rep_id)
          # load data
          if (covariate_name %in% c('missing_rate', 'pattern_gap', 'sequencing_depth')){covariate_val_format = format(round(covariate_val, 2),nsmall=1)} 
          else {covariate_val_format = covariate_val}
          
          file_name_zscore = paste0(folder, 'pairwise_zscores_combined_', covariate_key, '_', covariate_val_format, '_rep', rep_id, '.txt')
          file_name_meta = paste0(folder, 'pairwise_contrasts_metadata_', covariate_key, '_', covariate_val_format, '_rep', rep_id, '.txt')
          
          df_list = preprocess_input(file_name_zscore, file_name_meta)
          df_combined = df_list$df_combined
          df_meta = df_list$df_meta
          genes = rownames(df_combined)
          
          # run FClust 
          input = prepare_fpca(df_combined, df_meta, genes = genes, cell_type = cell_type, perturb = perturb)
          fpca_obj = FPCA(input$Ly, input$Lt)
          cout = FClust(input$Ly, input$Lt, k = 3, cmethod = cmethod,
                        optnsFPCA = fpca_obj$optns)
          ARI = calculate_cluster_accuracy(cout, n_patterns = 3)
          
          # run correlation measurement
          timepoints = seq(0,1,1/(n_timepoints-1))
          if (covariate_key %in% c('missing_rate', 'n_times', 'noise','sequencing_depth')){
            true_vals = retrieve_true_labels(type_pattern, n_timepoints, n_genes_pattern, mean_pert, mean_type_pert, pattern_coeff_gap)
            corr = calculate_correlation(timepoints = timepoints, fpca_obj = fpca_obj, true_vals = true_vals)
          }
          else if(covariate_key == 'pattern_gap'){
            true_vals = retrieve_true_labels(type_pattern, n_timepoints, n_genes_pattern, mean_pert, mean_type_pert, covariate_val)
            corr = calculate_correlation(timepoints = timepoints, fpca_obj = fpca_obj, true_vals = true_vals)
          } else { # 'nonlinear'
            corr = 0
          }
          
          # curate results
          result = c('', covariate_val, rep_id, cell_type, perturb, ARI, corr, 'FPCA_EMCluster')
          df_metrics[dim(df_metrics)[1]+1, ] = result
        }
      }
    }

  }
  return (df_metrics)
}

#### 
run_pipeline_over_DEGs = function(df_deg, perturb = 'Perturb_0', cell_type = 'Type_0', cmethod = 'EMCluster', 
                                  n_reps = 3, pattern_coeff_gap = 1.0, mean_pert = 0.3, 
                                  mean_type_pert = 2.0, n_patterns = 3, n_genes_pattern = 20, n_timepoints = 20)
{
  # initialization
  df_deg = df_deg[(df_deg$perturb == perturb) & (df_deg$cell_type == cell_type), ]
  n_genes = length(unique(df_deg$names))
  
  # add a column for gene ordering
  methods = c('t-test', 'wilcoxon')
  genes_order = 1:60
  genes = paste0('Gene_', 0:(n_genes-1))
  names(genes_order) = genes
  df_deg['genes_order'] = genes_order[df_deg$names]
  
  df_metrics = data.frame(matrix(ncol = 8, nrow = 0))
  colnames(df_metrics) = c('covariate_name', 'covariate_level', 'rep', 'cell_type', 'perturb', 'ARI', 'correlation', 'method')
  covariate_names = unique(df_deg$covariate_name)
  
  for (method in methods)
  {
    print(method)
    df_deg_sub1 = df_deg[df_deg[['method']] == method, ]
    for (covariate_name in covariate_names)
    {
      df_deg_sub2 = df_deg_sub1[df_deg_sub1[['covariate_name']] == covariate_name, ]
      covariate_values = unique(df_deg_sub2$covariate_level)
      for (covariate_val in covariate_values)
      {
        df_deg_sub3 = df_deg_sub2[df_deg_sub2[['covariate_level']] == covariate_val, ]
        if ((covariate_name == 'missing_rate') & (covariate_val == 0.7)) {next} # A stable solution is not avaliable for EMCluster
        for (rep in 0:(n_reps - 1))
        {
          if ((rep %in% c(1,2)) & (covariate_name == 'sequencing_depth')){next} 
          rep_name = paste0('rep',rep)
          df_deg_sub4 = df_deg_sub3[df_deg_sub3[['rep']] == rep_name, ]
          df_deg_sub4 = df_deg_sub4[with(df_deg_sub4, order(genes_order, time)), ]
          input = MakeFPCAInputs(IDs = df_deg_sub4$names,
                                 tVec = df_deg_sub4$time,
                                 yVec = df_deg_sub4$scores)
          fpca_obj = FPCA(input$Ly, input$Lt)
          cout = FClust(input$Ly, input$Lt, k = 3, cmethod = cmethod,
                        optnsFPCA = fpca_obj$optns)
          ARI = calculate_cluster_accuracy(cout, n_patterns = 3)
          
          # run correlation measurement
          timepoints = seq(0,1,1/(n_timepoints-1))
          if (covariate_name %in% c('missing_rate', 'n_times', 'noise','sequencing_depth')){
            true_vals = retrieve_true_labels(type_pattern, n_timepoints, n_genes_pattern, mean_pert, mean_type_pert, pattern_coeff_gap)
            corr = calculate_correlation(timepoints = timepoints, fpca_obj = fpca_obj, true_vals = true_vals)
          }
          else if(covariate_name == 'pattern_gap'){
            true_vals = retrieve_true_labels(type_pattern, n_timepoints, n_genes_pattern, mean_pert, mean_type_pert, covariate_val)
            corr = calculate_correlation(timepoints = timepoints, fpca_obj = fpca_obj, true_vals = true_vals)
          } else { # 'nonlinear'
            corr = 0
          }          
          # append results
          if (covariate_name %in% c('missing_rate', 'pattern_gap', 'sequencing_depth')){covariate_val_format = format(round(covariate_val, 2),nsmall=1)} 
          else {covariate_val_format = format(round(covariate_val, 2),nsmall=0)}
          
          result = c(covariate_name, covariate_val_format, rep, cell_type, perturb, ARI, corr, paste0('FPCA_EMCluster_', method))
          df_metrics[dim(df_metrics)[1]+1, ] = result
        }
      }
    }
  }
  return (df_metrics)
}

# run real benchmark datasets 
# 1. benchmark performance of fPCA across 6 covariates
folders = c('../simulation_missing_rate_',
            '../simulation_n_times',
            '../simulation_noise_level',
            '../simulation_pattern_gap_',
            '../simulation_sequencing_depth_',
            '../../simulation_nonLinear/simulation_nonlinear_covariates')
covariate_names = c('missing_rate', 'n_times', 'noise', 'pattern_gap', 'sequencing_depth', 'nonlinear_covariate_level')
covariate_values_vector = list('missing_rate' = seq(0,0.7,0.1),
                               'n_times' = seq(5,45,5),
                               'noise' = seq(0,9,1),
                               'pattern_gap' = seq(0.5,3.8,0.3),
                               'sequencing_depth' = seq(0.7,3.7,0.3),
                               'nonlinear_covariate_level' = seq(1,3,1))

df_metrics_all = data.frame(matrix(ncol = 8, nrow = 0))
colnames(df_metrics_all) = c('covariate_name', 'covariate_level', 'rep', 'ARI', 'cell_type', 'perturb', 'correlation', 'method')

for (i in 1:length(covariate_names))
{
  folder = folders[i]
  covariate_name = covariate_names[i]
  covariate_values = covariate_values_vector[[covariate_name]]
  print(covariate_name)
  
  df_metrics = run_pipeline_over_zscore(covariate_key = covariate_name, covariate_values = covariate_values,
                                        folder = paste0(folder, '/glm_zscore/'), pattern_coeff_gap = 1.0, mean_pert = 0.3, 
                                        mean_type_pert = 2.0, n_patterns = 3, n_genes_pattern = 20, n_timepoints = 20, n_rep = 3,
                                        cmethod = 'EMCluster')
  df_metrics['covariate_name'] = covariate_name
  df_metrics_all = rbind(df_metrics_all, df_metrics)
}
write.table(df_metrics_all, 'results_clusters_fPCA_EMCluster.txt', sep = '\t')


# 2. benchmark performance of DE scores across 6 covariates
df_DEG = read.table('../DEGs/DEG_concat.txt', sep = '\t', header = 1)
df_DEG$covariate_level = gsub('_','', df_DEG$covariate_level)
df_DEG$covariate_level = as.numeric(df_DEG$covariate_level)

df_metrics_all_deg = data.frame(matrix(ncol = 8, nrow = 0))
colnames(df_metrics_all_deg) = c('covariate_name', 'covariate_level', 'rep', 'ARI', 'cell_type', 'perturb', 'correlation', 'method')

for (cell_type in c('Type_0', 'Type_1'))
{
  for (perturb in c('Perturb_0', 'Perturb_1'))
  {
    df_metrics_all_deg_sub = run_pipeline_over_DEGs(df_DEG, perturb = perturb, cell_type = cell_type, cmethod = 'EMCluster', 
                                                n_reps = 3, pattern_coeff_gap = 1.0, mean_pert = 0.3, 
                                                mean_type_pert = 2.0, n_patterns = 3, n_genes_pattern = 20)
    df_metrics_all_deg = rbind(df_metrics_all_deg, df_metrics_all_deg_sub)
  }
}
write.table(df_metrics_all_deg, 'results_clusters_DEG_fPCA_EMCluster.txt', sep = '\t')



