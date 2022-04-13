library(stringr)

setwd('/Volumes/aronow/Kang/CellDrift/test_examples/covid_blood_atlas/glm_time_v3_all/fpca/')
df_combined = read.table('pairwise_zscores_HVG1000_combined_all.txt', sep = '\t', header = 1, row.names = 1)

## From Agresti(2007) p.39
get_qualified_genes = function(cell_type, perturb, rep = 'rep1')
{
  # get files
  path = paste0('../', rep, '/output_celldrift/')
  target_files = c()
  timepoints = c()
  for (file_name in list.files(path)){
    if (grepl('pairwise', file_name) & grepl(perturb, file_name)){
      target_files = append(target_files, file_name)
      timepoints = append(timepoints, as.integer(str_replace(strsplit(file_name, '_')[[1]][5], 'd','')))
    }
  }
  print(target_files)
  # get genes
  genes = c()
  new_type = paste0("(", cell_type, ", ", cell_type, ")")
  for (file_name in target_files){
    df = read.table(paste0(path, file_name), sep = '\t', header = 1)
    print(unique(df$cell_type))
    df = df[df$cell_type == new_type, ]
    sign_genes = df[df$p_fdr < 0.05, 'gene']
    genes = append(genes, sign_genes)
  }
  return(unique(genes))
}

test3 = get_qualified_genes(cell_type = 'cMono', perturb = 'COVID_CRIT', rep = 'rep3')

overlap_genes = intersect(intersect(test, test2), test3)
write.table(data.frame('names' = overlap_genes, 'dummy' = 1), 'selected_overlap_features.txt', sep = '\t')
saveRDS(overlap_genes, 'selected_overlap_features.rds')
