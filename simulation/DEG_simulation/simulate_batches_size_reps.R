library(Seurat)
library(splatter)
library(scater)
library(SeuratDisk)

ifn_data = readRDS("batch2_all_raw_filtered_v2.rds")

## simulate data for specific cell types
# extract parameters of CD14+ Monocytes
ifn_mono = ifn_data[,ifn_data$cell  == "CD14+ Monocytes"]
ifn_mono = as.SingleCellExperiment(ifn_mono)
param_mono = splatEstimate(as.matrix(ifn_mono@assays@data$counts))
saveRDS(param_mono, 'param_estimates_cd14mono_splatter.rds')

runSeurat = function(seu)
{
  seu@meta.data[["nCount"]] = colSums(seu@assays$originalexp@counts)
  seu@meta.data[["nGene"]] = colSums(seu@assays$originalexp@counts > 0)
  seu@assays$originalexp@meta.features$sparity = rowSums(seu@assays$originalexp@counts == 0 ) / dim(seu)[2]
  seu = NormalizeData(seu)
  seu = FindVariableFeatures(seu, nfeatures = 400)
  seu = ScaleData(seu)
  seu = RunPCA(seu, features = VariableFeatures(seu))
  seu = FindNeighbors(seu)
  seu = RunUMAP(seu, dim = 1:30)
  return(seu)
}

define_params_run_seurat = function(param_mono,
                                    seed = 42,
                                    nGenes = 3000,
                                    nGroup = 2, 
                                    group.prob = c(0.5, 0.5),
                                    de.prob = 0.5,
                                    de.downProb = 0.5,
                                    nCells = 1000,
                                    nBatches = 3,
                                    batchCells = c(600, 200, 200),
                                    batch.facLoc = 0.1,
                                    batch.facScale = 0.1,
                                    out.prob = 0,
                                    mean.rate = 15,
                                    lib.loc = 8.5)
{
  # define group and DEG parameters 
  param_mono@seed = seed
  param_mono@nGenes = nGenes
  param_mono@nGroups = nGroup
  param_mono@group.prob = group.prob
  param_mono@de.prob = de.prob
  param_mono@de.downProb = de.downProb
  param_mono@nCells = nCells
  param_mono@nBatches = nBatches
  param_mono@batchCells = batchCells
  param_mono@batch.facLoc = batch.facLoc
  param_mono@batch.facScale = batch.facScale
  param_mono@out.prob = out.prob
  param_mono@mean.rate = mean.rate
  param_mono@lib.loc = lib.loc
  
  # simulate the data
  sim_mono = splatSimulate(param_mono, method = "groups")
  
  # create seurat object and do downstream analysis
  seu = as.Seurat(sim_mono, data = NULL)
  # introduce imbalance in some batches
  remove_cells = sample(colnames(seu[, (seu$Batch == "Batch1") & (seu$Group == "Group1")]), 260)
  remain_cells = setdiff(colnames(seu), remove_cells)
  
  seu = seu[, remain_cells]
  seu = runSeurat(seu)
  return(seu)
}

# run with three batch parameters
create_test_data_rep = function(rep_time = 10)
{
  for (x in 1:rep_time)
  {
    random_seeds = c(1,2,3,4) + x
    vct_batch_facloc = c(0.02, 0.1, 0.4, 0.7)
    vct_batch_facscale = c(0.1, 0.1, 0.1, 0.1)
    seu_list = list()
    print(random_seeds)
    for (i in 1:length(vct_batch_facloc))
    {
      print(i)
      batch_facloc = vct_batch_facloc[i]
      batch_facscale = vct_batch_facscale[i]
      random_seed = random_seeds[i]
      seu = define_params_run_seurat(param_mono, 
                                     seed = random_seed, 
                                     batch.facLoc = batch_facloc, 
                                     batch.facScale = batch_facscale)
      seu_list[[i]] = seu
    }
    
    # transform data types
    dir.create(paste0("Rep", x))
    for (j in 1:length(vct_batch_facloc)){
      print(j)
      batch_facloc = vct_batch_facloc[j]
      seu_list[[j]]@assays$originalexp@data = seu_list[[j]]@assays$originalexp@counts

      file_name = paste0("Rep", x, "/sim_batchLoc_",batch_facloc,".h5Seurat")
      SaveH5Seurat(seu_list[[j]], file_name)
      Convert(file_name, dest = "h5ad")
      file_name = paste0("Rep", x, "/sim_batchLoc_",batch_facloc,".rds")
      saveRDS(seu_list[[j]], file_name)
    }
  }
}

setwd("test_batch_h5ad_v2/")
create_test_data_rep(10)
