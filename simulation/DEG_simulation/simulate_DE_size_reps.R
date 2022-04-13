# various de.prob
library(Seurat)
library(splatter)
library(scater)
library(SeuratDisk)

# load splatter parameters estimated from a reference data
param_mono = readRDS("param_estimates_cd14mono_splatter.rds")

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
                                    nCells = 600,
                                    nBatches = 1,
                                    batchCells = c(600),
                                    batch.facLoc = 0.02,
                                    batch.facScale = 0.02,
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
  seu = runSeurat(seu)
  return(seu)
}

# run with multiple DE probabilities
create_test_data_rep = function(rep_time = 1)
{
  for (x in 1:rep_time)
  {
    random_seeds = c(1,2,3,4) + x
    de.prob.list = c(0.05, 0.2, 0.5, 0.8)
    seu_list = list()
    print(random_seeds)
    for (i in 1:length(de.prob.list))
    {
      print(i)
      de.prob = de.prob.list[i]
      random_seed = random_seeds[i]
      seu = define_params_run_seurat(param_mono, 
                                     seed = random_seed, 
                                     de.prob = de.prob)
      seu_list[[i]] = seu
    }
    
    # transform data types
    dir.create(paste0("Rep", x))
    for (j in 1:length(de.prob.list)){
      print(j)
      de.prob = de.prob.list[j]
      seu_list[[j]]@assays$originalexp@data = seu_list[[j]]@assays$originalexp@counts
      
      file_name = paste0("Rep", x, "/sim_deProb_",de.prob,".h5Seurat")
      saveRDS(seu_list[[j]], paste0("Rep", x, "/sim_deProb_",de.prob,".rds"))
      SaveH5Seurat(seu_list[[j]], file_name)
      Convert(file_name, dest = "h5ad")
    }
  }
}

setwd("test_deProb_h5ad_v1/")
create_test_data_rep(10)
