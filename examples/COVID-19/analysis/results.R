library(stringr)
library(fdapace)
library(EMCluster)
library(mclust)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)

setwd('/Volumes/aronow/Kang/CellDrift/test_examples/covid_blood_atlas/glm_time_v3_all/results/')
source('function_.R')

# 0. visualization setting
f = function(pal) brewer.pal(brewer.pal.info[pal, 'maxcolors'], pal)
cols = c()
for (pl in c('RdBu', 'PiYG')){
  cols = c(cols, f(pl))
}
cols = setdiff(cols, c('#F7F7F7'))

# 1. load data
selected_genes = readRDS('../fpca/selected_overlap_features.rds')

cout = readRDS('../fpca/selected_genes_20_EMClust_cluster.rds')
fpca = cout$fpca

df_cluster = read.table('../fda_analysis/fda_clusters.txt', sep = '\t', header = 1, row.names = 1)
df_cluster$clusters_EMCluster = cout$cluster

df_combined = read.table('../fpca/pairwise_zscores_HVG1000_combined_all.txt', sep = '\t', header = 1, row.names = 1)
df_combined[is.na(df_combined)] = 0
df_meta = read.table('../factor_analysis/metadata_combined.txt', sep = '\t', header = 1, row.names = 1)
rownames(df_meta) = str_replace(rownames(df_meta), '-', '.')
df_meta$time_v2 = as.integer(str_replace(df_meta$time, 'd',''))
df_meta$id = paste0(df_meta$cell_type, '_', df_meta$perturbation)

# 2. define functions
visualize_cluster_fda = function(fpca_obj, cout, col, cluster = 1, k = 3)
{
  clusters = cout$cluster
  subset_ids = which(cout$cluster == cluster)
  p = CreatePathPlot(fpca_obj, subset = subset_ids, K = k, showMean = T, showObs = F,col = col,
                     main = paste0('Fuzzy-kmeans Cluster ', cluster, ' of classical Monocytes in Severe COVID-19'), 
                     lty = 1, pch = 4, xlab = 'days since onset', ylab = 'contrast coefficients'); grid()
  return (p)
}

# 3. basic visualization
pdf('figures/FPCA_basic_cMono_Sev.pdf', width = 10, height = 10)
plot(fpca)
dev.off()

# 4. draw clusters
cout$cluster = df_cluster$clusters_fuzzy
for (cluster in 0:19)
{
  col = cols[cluster + 1]
  pdf(paste0('figures/curves_cluster/cluster_curves_', cluster, '.pdf'), width = 5, height = 5)
  visualize_cluster_fda(fpca, cout, cluster = cluster, k = 3, col)
  dev.off()
}

# 5. time warping
query_genes = c('IFI27', 'IFITM2','IFITM3', 'IFI6')
query_genes = c('CTSD')
cell_types = c('cMono')
perturbations = c('COVID_HCW_MILD', 'COVID_MILD', 'COVID_SEV', 'COVID_CRIT', 'Flu', 'Sepsis')
colours = c('#CC0000', '#4C9900', '#00CCCC', '#FF8800', '#6666FF', '#7F00FF')
colours_v2 = rep(colours, each = 3)
run_draw_curves_allReps(df_combined, df_meta, query_genes, cell_types = cell_types,
                       perturbations, colours_v2, template_perturbation = 'COVID_SEV',
                       step_pattern = asymmetric, open_end = T)
for (rep in c('rep1', 'rep2', 'rep3'))
{
  run_draw_curves_perRep(df_combined, df_meta, query_genes, cell_types = cell_types,
                         perturbations, colours, template_perturbation = 'COVID_SEV',
                         step_pattern = asymmetric, open_end = T, rep = rep)
}

# 5.2 trimmed alignment
cell_types = c('cMono')
perturbations = c('COVID_HCW_MILD', 'COVID_MILD', 'COVID_SEV', 'COVID_CRIT', 'Flu', 'Sepsis')
# df_meta_trim = df_meta[(df_meta$time_v2 < 20) & (df_meta$time_v2 > 1), ]
df_meta_trim = df_meta[(df_meta$time_v2 < 20), ]
df_combined_trim = df_combined[, rownames(df_meta_trim)]

query_genes = c('S100A8', 'S100A9', 'S100A12', 'CTSD', 'LGALS1')
for (query_gene in query_genes)
{
  for (rep in c('rep1', 'rep2', 'rep3'))
  {
    # run_draw_curves_perRep(df_combined_trim, df_meta_trim, query_gene, cell_types = cell_types,
    #                        perturbations, colours, template_perturbation = 'COVID_SEV',
    #                        step_pattern = asymmetric, open_end = T, rep = rep, validation = F)
    run_draw_curves_perRep(df_combined_trim, df_meta_trim, query_gene, cell_types = cell_types,
                           perturbations, colours, template_perturbation = 'COVID_SEV',
                           step_pattern = asymmetric, open_end = F, rep = rep, validation = F)
  }
}

# 6. One way anova
# identify anova significant genes
cell_types = c('cMono')
perturbations = c('COVID_HCW_MILD', 'COVID_MILD', 'COVID_SEV', 'COVID_CRIT', 'Flu', 'Sepsis')
colours = c('#CC0000', '#4C9900', '#00CCCC', '#FF8800', '#6666FF', '#7F00FF')
df_meta_trim = df_meta[(df_meta$time_v2 < 20), ]
df_combined_trim = df_combined[, rownames(df_meta_trim)]

df_anova = read.table('anova_results/combined_anova_results.txt', sep = '\t', header = 1, row.names = 1)
df_anova = df_anova[order(-df_anova$Statistic), ]
up_genes = rownames(df_anova[df_anova$mean.gap..severe...mild. > 0, ])[1:5]
dn_genes = head(rownames(df_anova[df_anova$mean.gap..severe...mild. < 0,]), 5)
plain_genes = tail(rownames(df_anova), 5)
for (query_gene in plain_genes)
{
  for (rep in c('rep1', 'rep2', 'rep3'))
  {
    run_draw_curves_perRep(df_combined_trim, df_meta_trim, query_gene, cell_types = cell_types,
                           perturbations, colours, template_perturbation = 'COVID_SEV',
                           step_pattern = asymmetric, open_end = F, rep = rep, validation = F)
  }
}


# 7. cell-type comp

# 8. Validation 
# 8.1 load data
df_val_combined = read.table('/Volumes/aronow/Kang/CellDrift/test_examples/covid_atlas_Ren_validation/glm_local/fda_analysis/pairwise_zscores_combined.txt',
                             sep = '\t', header = 1, row.names = 1)
df_val_meta = read.table('/Volumes/aronow/Kang/CellDrift/test_examples/covid_atlas_Ren_validation/glm_local/fda_analysis/pairwise_contrasts_metadata.txt',
                         sep = '\t', header = 1, row.names = 1)
df_val_combined[is.na(df_val_combined)] = 0
rownames(df_val_meta) = str_replace(rownames(df_val_meta), '-', '.')
df_val_meta$time_v2 = as.integer(str_replace(df_val_meta$time, 'd',''))
df_val_meta$id = paste0(df_val_meta$cell_type, '_', df_val_meta$perturbation, '_', df_val_meta$rep)

# 8.2 draw curves
query_genes = c('CTSD')
cell_types = c('Mono')
perturbations = c('mild', 'severe')

colours = c('#4C9900','#CC0000')
colours_v2 = rep(colours, each = 3)
run_draw_curves_allReps(df_val_combined, df_val_meta, query_genes, cell_types = cell_types,
                        perturbations, colours_v2, template_perturbation = 'severe',
                        step_pattern = asymmetric, open_end = T, validation = T)

# 8.3 trim timepoints
query_genes = c('S100A8', 'S100A9', 'S100A12', 'CTSD', 'LGALS1')
df_val_meta_trim = df_val_meta[df_val_meta$time_v2 < 20, ]
df_val_combined_trim = df_val_combined[, rownames(df_val_meta_trim)]
for (query_gene in query_genes)
{
  run_draw_curves_allReps(df_val_combined_trim, df_val_meta_trim, query_gene, cell_types = cell_types,
                          perturbations, colours_v2, template_perturbation = 'severe',
                          step_pattern = asymmetric, open_end = T, validation = T)
}



