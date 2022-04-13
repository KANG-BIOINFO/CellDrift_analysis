library(stringr)
library(fdapace)
library(EMCluster)
library(dtw)
library(ggplot2)
library(cowplot)

setwd('../test_examples/covid_blood_atlas/glm_time_v3_all/fpca/')

# prepare fpca input
prepare_fpca = function(df_zscore, df_meta, genes, cell_types, perturbs, rep = 'rep1'){
  # prepare input data 
  df_submeta = df_meta[(df_meta$cell_type %in% cell_types) & (df_meta$perturbation %in% perturbs), ]
  rep_specific = rownames(df_submeta)[grepl(rownames(df_submeta), pattern = rep)]
  df_submeta = df_submeta[rep_specific, ]
  df_submeta = df_submeta[order(df_submeta$time_v2),]

  comparisons = rownames(df_submeta)
  df_zscore = df_zscore[genes, ]
  df_zscore = df_zscore[, comparisons]

  yVec = c(); tVec = c(); IDs = c()
  for (i in 1:ncol(df_zscore)){
    col = colnames(df_zscore)[i]
    yVec = append(yVec, df_zscore[,i])
    tVec = append(tVec, rep(df_submeta[col, 'time_v2'], nrow(df_zscore)))
    IDs = append(IDs, paste0(df_submeta[col, 'id'], '-', rownames(df_zscore)))
  }
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
      # lines(ly[alignment$index1]~alignment$index2, col = 'blue')
      
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
                          perturbations, colours, template_perturbation = 'COVID_SEV',
                          step_pattern = asymmetric, open_end = T, rep = 'rep1'){
  
  # initialize
  template_name = paste0(cell_types, '_', template_perturbation, '-', query_genes)
  
  # prepare input and dtw aligned samples
  input = prepare_fpca_allReps(df_combined, df_meta, query_genes, 
                       cell_types = cell_types, perturbs = perturbations, rep = rep)
  aligned = perform_alignment(input, template_name = template_name, step_patten = step_pattern,
                              open_end = open_end)
  
  # adjust index
  names_raw = gsub(paste0(cell_types, '_'), '', gsub(paste0('-', query_genes), '', unlist(input$Lid)))
  names_aligned = gsub(paste0(cell_types, '_'), '', gsub(paste0('-', query_genes), '', unlist(aligned$Lid)))

  subset_idx_raw = match(perturbations, names_raw)
  subset_idx_aligned = match(perturbations, names_aligned)
  
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

####################
# run on real data #
####################

# 1. load data
df_combined = read.table('pairwise_zscores_HVG1000_combined_all.txt', sep = '\t', header = 1, row.names = 1)
df_meta = read.table('../factor_analysis/metadata_combined.txt', sep = '\t', header = 1, row.names = 1)
rownames(df_meta) = str_replace(rownames(df_meta), '-', '.')
df_meta$time_v2 = as.integer(str_replace(df_meta$time, 'd',''))
df_meta$id = paste0(df_meta$cell_type, '_', df_meta$perturbation)

# 2. run 
query_genes = c('IFIT3', 'IFI16', 'IFI35', 'IFI44', 'IFI44L', 'IFIT1', 'IFIT2')
query_genes = c('RPS11', 'RPS12', 'RPS13')
query_genes = c('RETN', 'NFE2', 'HRH2', 'PKM')
query_genes = c('S100A12', 'S100A8', 'CD74', 'S100A9', 'LGALS2', 'TMEM176B', 
                'PLBD1', 'MTRNR2L8', 'APOBEC3A', 'SGK1', 'TMEM176A', 'NKG7', 
                'CTSD', 'GIMAP7', '')
# query_genes = c('S100A8', 'S100A9')
cell_types = c('cMono')
perturbations = c('COVID_HCW_MILD', 'COVID_MILD', 'COVID_SEV', 'COVID_CRIT', 'Flu', 'Sepsis')
colours = c('#4C9900', '#00CCCC', '#FF8800', '#CC0000', '#6666FF', '#7F00FF')

for (query_gene in query_genes)
{
  run_draw_curves(df_combined, df_meta, query_gene, cell_types, 
                  perturbations, colours, template_perturbation = 'COVID_SEV',
                  open_end = T)
}


