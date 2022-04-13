library(stringr)
library(fdapace)
library(EMCluster)
library(dtw)
library(ggplot2)
library(cowplot)

setwd('../test_examples/covid_atlas_Ren_validation/glm_local/fda_analysis/')

# prepare fpca input
prepare_fpca = function(df_zscore, df_meta, genes, cell_types, perturbs){
  # prepare input data 
  df_submeta_all = df_meta[(df_meta$cell_type %in% cell_types) & (df_meta$perturbation %in% perturbs), ]
  # collect information from replicates
  reps = unique(df_meta$rep)
  yVec = c(); tVec = c(); IDs = c()
  
  for (rep in reps)
  {
    # prep ( filtering )
    rep_specific = rownames(df_submeta_all)[grepl(rownames(df_submeta_all), pattern = rep)]
    df_submeta = df_submeta_all[rep_specific, ]
    df_submeta = df_submeta[order(df_submeta$time_v2),]
    
    # get comparisons and values
    comparisons = rownames(df_submeta)
    df_zscore_rep = df_zscore[genes, ]
    df_zscore_rep = df_zscore_rep[, comparisons]
    
    # append values
    for (i in 1:ncol(df_zscore_rep)){
      col = colnames(df_zscore_rep)[i]
      yVec = append(yVec, df_zscore_rep[,i])
      tVec = append(tVec, rep(df_submeta[col, 'time_v2'], nrow(df_zscore_rep)))
      IDs = append(IDs, paste0(df_submeta[col, 'id'], '-', rep, '-', rownames(df_zscore_rep)))
    }
  }
  # initialize a FPCA input object
  input = MakeFPCAInputs(IDs = IDs,
                         tVec = tVec,
                         yVec = yVec)
  return (input)
}

# perform alignment
perform_alignment = function(input, template_name, step_patten = asymmetric, open_end = T)
{
  # initialization
  template = input$Ly[[template_name]]
  template_lt = input$Lt[[template_name]]
  xs_aligned = c()
  ys_aligned = c()
  ID_aligned = c()
  
  xs_aligned = append(xs_aligned, template_lt)
  ys_aligned = append(ys_aligned, template)
  ID_aligned = append(ID_aligned, rep(template_name, 
                                      length(input$Lt[[template_name]])))
  
  # perform dtw
  for (i in 1:length(input$Lid)){
    id = input$Lid[[i]]
    if (id != template_name)
    { 
      lt = input$Lt[i][[1]]
      ly = input$Ly[i][[1]]
      
      alignment = dtw(x = ly, y = template, keep.internals = F, 
                      window.type = 'none', step.pattern = step_patten,
                      open.begin = F, open.end = open_end)
      
      # plot(template)
      new_ly = ly[alignment$index1]
      new_lt = template_lt[alignment$index2]
      xs_aligned = append(xs_aligned, new_lt)
      ys_aligned = append(ys_aligned, new_ly)
      ID_aligned = append(ID_aligned, rep(id, length(new_lt)))
    }
  }
  # make aligned FPCA input data
  aligned_material = list('IDs' = ID_aligned, 'tVec' = xs_aligned, 'yVec' = ys_aligned)
  aligned = MakeFPCAInputs(IDs = ID_aligned,
                           tVec = xs_aligned,
                           yVec = ys_aligned,
                           deduplicate = T)
  # return (aligned)
  return (aligned)
}

run_draw_curves = function(df_combined, df_meta, query_genes, cell_types,
                           perturbations, colours, template_perturbation = 'severe',
                           step_pattern = asymmetric, open_end = T, rep = 'rep1'){
  
  # initialize
  template_name = paste0(cell_types, '_', template_perturbation, '-rep0-', query_genes)
  
  # prepare input and dtw aligned samples
  input = prepare_fpca(df_combined, df_meta, query_genes, 
                       cell_types = cell_types, perturbs = perturbations)
  aligned = perform_alignment(input, template_name = template_name, step_patten = step_pattern,
                              open_end = open_end)
  
  # adjust index
  subset_idx_raw = match(sort(unlist(input$Lid)), unlist(input$Lid))
  subset_idx_aligned = match(sort(unlist(input$Lid)), unlist(aligned$Lid))
  
  # run fPCA
  fpca_raw = FPCA(input$Ly, input$Lt)
  fpca_aligned = FPCA(aligned$Ly, aligned$Lt)
  fpca_aligned_smallBwCov = FPCA(aligned$Ly, aligned$Lt, list('userBwCov' = 0.8))
  
  # draw plots
  pdf(paste0('figures/aligned_examples/FPCA_raw_', cell_types, '_', query_genes, '.pdf'))
  CreatePathPlot(fpca_raw, subset = subset_idx_raw, K = fpca_raw$selectK, 
                 col = colours, pch = '.', 
                 cex = 6, lty = 'solid', lwd = 1.5,
                 xlab = 'days since onset', ylab = 'contrast coefficient'); grid()
  legend(x = 18, y = max(unlist(input$Ly)), legend = perturbations, col = colours, lty = 1, cex = 0.85)
  dev.off()
  pdf(paste0('figures/aligned_examples/FPCA_aligned_', cell_types, '_', query_genes, '.pdf'))
  CreatePathPlot(fpca_aligned, subset = subset_idx_aligned, K = fpca_aligned$selectK,
                 col = colours, pch = '.', 
                 cex = 6, lty = 'solid', lwd = 1.5,
                 xlab = 'days since onset', ylab = 'contrast coefficient'); grid()
  legend(x = 18, y = max(unlist(aligned$Ly)), legend = perturbations, col = colours, lty = 1, cex = 0.85)
  dev.off()
  pdf(paste0('figures/aligned_examples/FPCA_aligned_smallCovBw_', cell_types, '_', query_genes, '.pdf'))
  CreatePathPlot(fpca_aligned_smallBwCov, subset = subset_idx_aligned, K = fpca_aligned_smallBwCov$selectK,
                 col = colours, pch = '.', 
                 cex = 6, lty = 'solid', lwd = 1.5,
                 xlab = 'days since onset', ylab = 'contrast coefficient'); grid()
  legend(x = 18, y = max(unlist(aligned$Ly)), legend = perturbations, col = colours, lty = 1, cex = 0.85)
  dev.off()
}

# preparation of files
df_zscores = read.table('pairwise_zscores_combined.txt', sep = '\t', header = 1, row.names = 1)
df_zscores[is.na(df_zscores)] = 0

df_meta = read.table('pairwise_contrasts_metadata.txt', sep = '\t', header = 1, row.names = 1)
rownames(df_meta) = str_replace(rownames(df_meta), '-', '.')
df_meta$time_v2 = as.integer(str_replace(df_meta$time, 'd',''))
df_meta$id = paste0(df_meta$cell_type, '_', df_meta$perturbation)

# run FPCA and alignment
query_genes = c('S100A8', 'S100A9', 'S100A12', 'CTSD', 'LGALS1')
# query_genes = c('S100A8', 'S100A9')
cell_types = c('Mono')
perturbations = c('mild', 'severe')
colours = rep(c('#4C9900', '#FF8800'), each = 3)

for (query_gene in query_genes)
{
  run_draw_curves(df_zscores, df_meta, query_gene, cell_types, 
                  perturbations, colours, template_perturbation = 'severe',
                  open_end = T)
}