library(Seurat)
library(SingleCellExperiment)
library(MAST)

setwd("/Volumes/aronow/Kang/CellDrift/test/test_benchmark/simulated_ifn_data/test_deProb_celldrift_others/")

# test
# param_mono = readRDS("../../benchmark_neighborhoods/param_estimates_cd14mono_splatter.rds")
# test_data = define_params_run_seurat(param_mono)
# test_sce = as.SingleCellExperiment(test_data)
# test_sca = SceToSingleCellAssay(test_sce)
# 
# zlmCond = zlm(~ Group + Batch, sca = test_sca)
# summaryCond <- summary(zlmCond, doLRT='GroupGroup2')
# 
# summaryDt <- summaryCond$datatable
# fcHurdle <- merge(summaryDt[contrast=='GroupGroup2' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
#                   summaryDt[contrast=='GroupGroup2' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
# 
# output = as.data.frame(fcHurdle)


source_path = "../test/test_benchmark/simulated_ifn_data/test_deProb_h5ad_v1/"
rep_dirs = setdiff(list.dirs(source_path), source_path)

for (rep_dir in rep_dirs)
{
  print(rep_dir)
  rep_name = strsplit(rep_dir, "//")[[1]][3]
  rds_files = list.files(rep_dir, pattern = ".rds")

  output_all = data.frame()
  for (rds_file in rds_files)
  {
    DE_facLoc = strsplit(rds_file, "_")[[1]][3]
    DE_facLoc = as.numeric(strsplit(DE_facLoc, ".rds")[[1]][1])
    
    seu = readRDS(paste0(rep_dir, "/", rds_file))
    seu = NormalizeData(seu)
    sim = as.SingleCellExperiment(seu)
    print(sim@assays@data$logcounts)
    sca = SceToSingleCellAssay(sim)
    
    print("Start to run hurdle model")
    zlmCond = zlm(~ Group, sca = sca) # no need to adjust batch effect here, since it's neglectable.
    summaryCond <- summary(zlmCond, doLRT='GroupGroup2')

    summaryDt <- summaryCond$datatable
    fcHurdle <- merge(summaryDt[contrast=='GroupGroup2' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                      summaryDt[contrast=='GroupGroup2' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
    output = as.data.frame(fcHurdle)
    output[["DE_FacLoc"]] = DE_facLoc
    output_all = rbind(output_all, output)
  }
  write.table(output_all, paste0("mast_results/", rep_name, "_mast_Group2_vs_Group1.txt"), sep = "\t")
}




