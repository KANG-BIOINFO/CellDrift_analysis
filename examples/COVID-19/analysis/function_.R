library(stringr)
library(fdapace)
library(EMCluster)
library(dtw)
library(ggplot2)
library(cowplot)

# prepare fpca input
prepare_fpca = function(df_zscore, df_meta, genes, cell_types, perturbs, rep = 'rep1'){
  # prepare input data
  df_submeta = df_meta[(df_meta$cell_type %in% cell_types) & (df_meta$perturbation %in% perturbs), ]
  rep_specific = rownames(df_submeta)[grepl(rownames(df_submeta), pattern = rep)]
  df_submeta = df_submeta[rep_specific, ]
  df_submeta = df_submeta[order(df_submeta$time_v2),]
  print(dim(df_submeta))

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

prepare_fpca_allReps = function(df_zscore, df_meta, genes, cell_types, perturbs){
  # prepare input data 
  df_submeta = df_meta[(df_meta$cell_type %in% cell_types) & (df_meta$perturbation %in% perturbs), ]
  df_submeta = df_submeta[order(df_submeta$time_v2),]
  print(dim(df_submeta))

  comparisons = rownames(df_submeta)
  df_zscore = df_zscore[genes, ]
  df_zscore = df_zscore[, comparisons]
  
  yVec = c(); tVec = c(); IDs = c()
  for (i in 1:ncol(df_zscore)){
    col = colnames(df_zscore)[i]
    yVec = append(yVec, df_zscore[,i])
    tVec = append(tVec, rep(df_submeta[col, 'time_v2'], nrow(df_zscore)))
    IDs = append(IDs, paste0(df_submeta[col, 'id'], '-' , rownames(df_zscore)))
    # IDs = append(IDs, paste0(df_submeta[col, 'id'], '-', tail(strsplit(col, '_')[[1]], 1),
    #                          '-' , rownames(df_zscore)))
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
  print(template_name)
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

visualize_genes_PCs = function(fpca_obj, genes, genes_all, pcs = c(1,2), fontsize = 0.5)
{
  df_pca = data.frame(fpca_obj$xiEst[, pcs])
  colnames(df_pca) = c('Component1', 'Component2')
  rownames(df_pca) = genes_all
  df_pca[['highlight']] = '0'
  df_pca[genes, 'highlight'] = '1'
  df_pca[['label']] = rownames(df_pca)
  
  p = ggplot(df_pca, aes(x = Component1, y = Component2, label=ifelse(highlight == '1', as.character(label),''))) +
    geom_point(aes(colour = highlight)) + 
    geom_text(size = fontsize) + 
    geom_text_repel(max.overlaps = 20) + geom_label_repel()+ 
    xlab(paste0('PC_', pcs[1])) + ylab(paste0('PC_', pcs[2]))  +
    theme(legend.position = 'none',
          panel.background = element_rect(fill = NA, colour = 'black'),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 12)) + 
    scale_colour_manual(values = c('0' = '#EEEEEE', '1' = '#CC0000'))
  
  return(p)
}

run_draw_curves_allReps = function(df_combined, df_meta, query_genes, cell_types,
                           perturbations, colours, template_perturbation = 'COVID_SEV',
                           step_pattern = asymmetric, open_end = T, validation = F){
  
  # initialize
  template_name = paste0(cell_types, '_', template_perturbation, '_rep0-', query_genes)
  print(template_name)
  # prepare input and dtw aligned samples
  input = prepare_fpca_allReps(df_combined, df_meta, query_genes, 
                              cell_types = cell_types, perturbs = perturbations)
  print(unlist(input$Lid))
  aligned = perform_alignment(input, template_name = template_name, step_patten = step_pattern,
                              open_end = open_end)
  
  # adjust index
  subset_idx_raw = match(sort(unlist(input$Lid)), unlist(input$Lid))
  subset_idx_aligned = match(sort(unlist(input$Lid)), unlist(aligned$Lid))
  
  print(sort(unlist(input$Lid)))
  print(subset_idx_raw)
  print(subset_idx_aligned)
  
  # run fPCA
  fpca_raw = FPCA(input$Ly, input$Lt)
  fpca_aligned = FPCA(aligned$Ly, aligned$Lt)
  # fpca_aligned_smallBwCov = FPCA(aligned$Ly, aligned$Lt, list('userBwCov' = 0.8))
  
  # draw plots
  name_folder = if (validation) 'figures/alignment_validation_trimmed/FPCA_raw_' else 'figures/alignment/FPCA_raw_'
  pdf(paste0(name_folder, cell_types, '_', query_genes, '.pdf'))
  CreatePathPlot(fpca_raw, subset = subset_idx_raw, K = fpca_raw$selectK, 
                 pch = '.', col = colours,
                 cex = 6, lty = 'solid', lwd = 2.5,
                 xlab = 'days since onset', ylab = 'contrast coefficient'); grid()
  legend(x = 15, y = max(unlist(input$Ly)), legend = sort(unlist(input$Lid)), 
         lty = 1, cex = 0.85, col = colours)
  dev.off()
  
  name_folder = if (validation) 'figures/alignment_validation_trimmed/FPCA_aligned_' else 'figures/alignment/FPCA_aligned_'
  pdf(paste0(name_folder, cell_types, '_', query_genes, '.pdf'))
  CreatePathPlot(fpca_aligned, subset = subset_idx_aligned, K = fpca_aligned$selectK,
                 col = colours, pch = '.', 
                 cex = 6, lty = 'solid', lwd = 2.5,
                 xlab = 'days since onset', ylab = 'contrast coefficient'); grid()
  legend(x = 15, y = max(unlist(aligned$Ly)), legend = sort(unlist(input$Lid)),
         col = colours, lty = 1, cex = 0.85)
  dev.off()
}

run_draw_curves_perRep = function(df_combined, df_meta, query_genes, cell_types,
                                   perturbations, colours, template_perturbation = 'COVID_SEV',
                                   step_pattern = asymmetric, open_end = T, rep = 'rep1', validation = F){
  
  # initialize
  template_name = paste0(cell_types, '_', template_perturbation, '-', query_genes)
  print(template_name)
  
  # prepare input and dtw aligned samples
  input = prepare_fpca(df_combined, df_meta, query_genes, 
                       cell_types = cell_types, perturbs = perturbations, rep = rep)
  print(unlist(input$Lid)[1:5])
  print(length(input$Lid))
  aligned = perform_alignment(input, template_name = template_name, step_patten = step_pattern,
                              open_end = open_end)
  
  # adjust index
  subset_idx_raw = match(sort(unlist(input$Lid)), unlist(input$Lid))
  subset_idx_aligned = match(sort(unlist(input$Lid)), unlist(aligned$Lid))
  
  print(sort(unlist(input$Lid)))
  print(subset_idx_raw)
  print(subset_idx_aligned)
  
  # run fPCA
  fpca_raw = FPCA(input$Ly, input$Lt)
  fpca_aligned = FPCA(aligned$Ly, aligned$Lt)
  # fpca_aligned_smallBwCov = FPCA(aligned$Ly, aligned$Lt, list('userBwCov' = 0.8))
  
  # draw plots
  name_folder = if (validation) 'figures/alignment_validation/FPCA_raw_' else 'figures/alignment_trimmed/FPCA_raw_'
  pdf(paste0(name_folder, cell_types, '_', query_genes, '_', rep, '.pdf'))
  CreatePathPlot(fpca_raw, subset = subset_idx_raw, K = fpca_raw$selectK, 
                 pch = '.', col = colours,
                 cex = 6, lty = 'solid', lwd = 3.5,
                 xlab = 'days since onset', ylab = 'contrast coefficient'); grid()
  legend(x = 15, y = max(unlist(input$Ly)), legend = sort(unlist(input$Lid)), 
         lty = 1, cex = 0.85, col = colours)
  dev.off()
  name_folder = if (validation) 'figures/alignment_validation/FPCA_aligned_' else 'figures/alignment_trimmed/FPCA_aligned_'
  pdf(paste0(name_folder, cell_types, '_', query_genes, '_', rep, '.pdf'))
  CreatePathPlot(fpca_aligned, subset = subset_idx_aligned, K = fpca_aligned$selectK,
                 col = colours, pch = '.', 
                 cex = 6, lty = 'solid', lwd = 3.5,
                 xlab = 'days since onset', ylab = 'contrast coefficient'); grid()
  legend(x = 15, y = max(unlist(aligned$Ly)), legend = sort(unlist(input$Lid)),
         col = colours, lty = 1, cex = 0.85)
  dev.off()
}