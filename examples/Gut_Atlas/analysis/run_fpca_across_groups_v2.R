library(stringr)
library(fdapace)
library(EMCluster)
library(dtw)
library(ggplot2)
library(cowplot)
library(RColorBrewer)

setwd('../test_examples/gut_development/glm_fda/glm_zscore/')

# prepare file for functional pca (subjects are genes)
prepare_fpca = function(df_zscore, df_meta, genes, cell_type = 'colon', perturb = 'mesenchymal'){
  # prepare input data 
  df_submeta = df_meta[(df_meta$cell_type == cell_type) & (df_meta$perturbation == perturb), ]
  df_submeta = df_submeta[order(df_submeta$time_v2),]
  print(dim(df_submeta))
  
  comparisons = rownames(df_submeta)
  df_zscore = df_zscore[genes, ]
  df_zscore = df_zscore[, comparisons]
  print(dim(df_submeta))
  
  yVec = c(); tVec = c(); IDs = c()
  for (i in 1:ncol(df_zscore)){
    col = colnames(df_zscore)[i]
    yVec = append(yVec, df_zscore[,i])
    tVec = append(tVec, rep(df_submeta[col, 'time_v2'], nrow(df_zscore)))
    IDs = append(IDs, rownames(df_zscore))
  }
  print(head(IDs))
  print(head(tVec))
  print(head(yVec))
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


# visualize cluster
visualize_cluster_fda = function(fpca_obj, cout, colour, cluster = 1, k = 3)
{
  clusters = cout$cluster
  subset_ids = which(cout$cluster == cluster)
  p = CreatePathPlot(fpca_obj, subset = subset_ids, K = k, 
                     showMean = T, showObs = F, 
                     main = paste0('Cluster_', cluster), 
                     lty = 'solid', col = colour, 
                     pch = 4, xlab = 'days since onset', ylab = 'contrast coefficients'); grid()
  return (p)
}

# visualize genes
visualize_genes_fda = function(fpca_obj, genes_all, genes_query, k = 4)
{
  genes_ids = c()
  for (gene in genes_query){genes_ids = append(genes_ids, which(genes_all == gene))}
  CreatePathPlot(fpca_obj, subset = genes_ids, K = k, main = paste0('K = ', k) , pch = 4); grid()
}

# visualize genes on principal components
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

get_clusters = function(cout, query_genes, all_genes){
  df_cluster = data.frame(cout$cluster)
  rownames(df_cluster) = all_genes
  colnames(df_cluster) = c('cluster')
  df_cluster = df_cluster[query_genes,]
  return(df_cluster)
}

run_draw_curves = function(df_combined, df_meta, query_genes, cell_types,
                          perturbations, colours, template_perturbation = 'COVID_SEV',
                          step_pattern = asymmetric, open_end = T, rep = 'rep1'){
  
  # initialize
  template_name = paste0(cell_types, '_', template_perturbation, '-', query_genes)
  
  # prepare input and dtw aligned samples
  input = prepare_fpca(df_combined, df_meta, query_genes, 
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

# 0. visualization setting
f = function(pal) brewer.pal(brewer.pal.info[pal, 'maxcolors'], pal)
cols = c()
for (pl in c('RdBu', 'PiYG')){
  cols = c(cols, f(pl))
}
cols = setdiff(cols, c('#F7F7F7'))

# 1. load data
df_combined = read.table('pairwise_zscores_combined_gut_development.txt', sep = '\t', header = 1, row.names = 1)
df_combined[is.na(df_combined)] = 0
df_meta = read.table('pairwise_contrasts_metadata_gut_development.txt', sep = '\t', header = 1, row.names = 1)
rownames(df_meta) = str_replace(rownames(df_meta), '-', '.')
df_meta$time_v2 = as.double(str_replace(df_meta$time, 'F',''))
df_meta$id = paste0(df_meta$cell_type, '_', df_meta$perturbation)

# 2. EMClusters
df_genes = read.table('significant_genes.txt', sep = '\t', header = 1, row.names = 1)
# cell_types = unique(df_meta$cell_type)
cell_types = c('epithelium', 'mesenchymal')
perturbations = unique(df_meta$perturbation)

n_clusters = 20
for (cell_type in cell_types)
{
  print(cell_type)
  for (perturbation in perturbations)
  {
    print(perturbation)
    # initialization
    group = paste0('(', cell_type, ', ', perturbation, ')')
    print(group)
    genes_selected = df_genes[df_genes$group == group, 'genes']
    
    # run functional PCA
    input = prepare_fpca(df_combined, df_meta, genes = genes_selected, cell_type = cell_type, perturb = perturbation)
    fpca_obj = FPCA(input$Ly, input$Lt)
    
    # find temporal clusters
    cout = FClust(input$Ly, input$Lt, k = n_clusters,
                  optnsFPCA = fpca_obj$optns)
    
    # output figures
    dir.create(paste0('figures/', cell_type, '_', perturbation))
    for (cluster_id in 1:20)
    {
      pdf(paste0('figures/', cell_type, '_', perturbation, '/EMCluster_', cluster_id, '.pdf'))
      visualize_cluster_fda(fpca_obj, cout, cluster = cluster_id, k = fpca_obj$selectK)
      dev.off()
    }
    
    # output clusters
    df_cout = data.frame('names' = unlist(input$Lid), 'clusters' = cout$cluster)
    write.table(df_cout, paste0('selected_genes_', n_clusters, '_EMClust_clusters_', cell_type, '_', perturbation, '.txt'), sep = '\t')
    saveRDS(cout, paste0('selected_genes_', n_clusters, '_EMClust_cluster_', cell_type, '_', perturbation, '.rds'))
  }
}
save.image('gut_glm.RData')

# 2.1 polish Clusters
n_clusters = 20
files_cout = c('selected_genes_20_EMClust_cluster_epithelium_colon.rds',
               'selected_genes_20_EMClust_cluster_epithelium_ileum.rds',
               'selected_genes_20_EMClust_cluster_mesenchymal_colon.rds',
               'selected_genes_20_EMClust_cluster_mesenchymal_ileum.rds')
folder_names = c('epithelium_colon', 'epithelium_ileum',
                 'mesenchymal_colon', 'mesenchymal_ileum')
for (i in 1:4)
{
  cout = readRDS(files_cout[i])
  folder_name = folder_names[i]
  fpca_obj = cout$fpca
  
  dir.create(paste0('figures/EMClusters_v2/', folder_name))
  
  for (cluster in 1:n_clusters)
  {
    colour = cols[cluster]
    pdf(paste0('figures/EMClusters_v2/', folder_name, '/EMCluster_', cluster, '.pdf'),
        width = 6, height = 6)
    visualize_cluster_fda(fpca_obj, cout, colour = colour, cluster = cluster, k = fpca_obj$selectK)
    dev.off()
  }
}


# 3. fuzzy cluster visualization
files_fuzzyclus = c('clusters_epithelium_colon.txt',
                    'clusters_epithelium_ileum.txt',
                    'clusters_mesenchymal_colon.txt',
                    'clusters_mesenchymal_ileum.txt')
cout_files = c('selected_genes_20_EMClust_cluster_epithelium_colon.rds',
               'selected_genes_20_EMClust_cluster_epithelium_ileum.rds',
               'selected_genes_20_EMClust_cluster_mesenchymal_colon.rds',
               'selected_genes_20_EMClust_cluster_mesenchymal_ileum.rds')
cell_types = c('epithelium','epithelium','mesenchymal','mesenchymal')
tissues = c('colon','ileum','colon','ileum')

n_clusters = 20
for (i in 1:4)
{
  df_fc = read.table(paste0('fuzzy_clusters/',files_fuzzyclus[i]), sep = '\t', header = 1, row.names = 1)
  cout = readRDS(cout_files[i])
  cell_type = cell_types[i]
  tissue = tissues[i]
  
  fc = df_fc$clusters_fuzzy
  cout$cluster = fc
  fpca_obj = cout$fpca
  
  # output figures
  print('create folder')
  dir.create(paste0('figures/fuzzy_clusters/', cell_type, '_', tissue))
  for (cluster_id in 0:(n_clusters-1))
  {
    if (cluster_id %in% unique(fc))
    {
      pdf(paste0('figures/fuzzy_clusters/', cell_type, '_', tissue, '/FuzzyCluster_', cluster_id, '.pdf'))
      visualize_cluster_fda(fpca_obj, cout, cluster = cluster_id, k = fpca_obj$selectK)
      dev.off()
    }
  }
}


# 4. show peaks of FPCA
# 4.1 epithelium (colon vs. duojejunum)
fpca = cout1$fpca
pdf('figures/fpca_basis_colon_epithelium.pdf', width = 6, height = 6)
plot(fpca)
dev.off()

# 4.3 mesenchymal (colon vs. duojejunum)
fpca3 = cout3$fpca
pdf('figures/fpca_basis_colon_mesenchymal.pdf', width = 6, height = 6)
plot(fpca3)
dev.off()
features_3 = read.table('selected_genes_20_EMClust_clusters_mesenchymal_colon.txt', 
                        sep = '\t', header = 1)[, 'names']
xiEst3 = data.frame(fpca3$xiEst)
rownames(xiEst3) = features_3
colnames(xiEst3) = paste0('FPC_', 1:2)
xiEst3$cluster = cout3$cluster
write.table(xiEst3, 'results/table_fpca_clusters_mesenchymal_colon.txt', sep = '\t')

save.image('gut_glm_v2.RData')


# plot FPCA-basic
pdf(file = 'figures/fpca_basis_ileum_epithelium.pdf')
plot(cout2$fpca)
dev.off()

pdf(file = 'figures/fpca_basis_ileum_mesenchymal.pdf')
plot(cout4$fpca)
dev.off()


