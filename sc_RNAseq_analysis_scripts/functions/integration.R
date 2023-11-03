library("Seurat")

# IntegratedSeuratAnalysis ----
## Standard 'Seurat' integration and analysis of cells in a list of Seurat objects

IntegratedSeuratAnalysis <- function(seurat_filtered,
                                     n_variable_feats,
                                     npcs_integration,
                                     npcs_clustering) {
  
  require(Seurat)
  require(scran)
  require(dplyr)
  
  
  if (!is.list(seurat_filtered)) {
    
    seurat_filtered <- Seurat::SplitObject(seurat_filtered, split.by = "orig.ident")
    
  }
  
  
  select <- dplyr::select
  
  #load("C:/Users/Fotini/Documents/Meduoa_research/scrna_seq/data/cycle.rda")#cell cycle genes must be loaded from a file
  
  
  for (seurat in 1:length(seurat_filtered)) {
    
    tryCatch({
      
      set.seed(1993)
      
      message("Analysis of sample: ")
      
      print(seurat_filtered[[seurat]]$orig.ident[[1]])
      
      message("Normalizing...")
      
      seurat_filtered[[seurat]] <- Seurat::NormalizeData(seurat_filtered[[seurat]])#, assay = NULL, normalization.method = "LogNormalize", scale.factor = 10000, margin = 1, verbose = F)
      
      #message("Cell Cycle Scoring...")
      
      #seurat_filtered[[seurat]] <- Seurat::CellCycleScoring(seurat_filtered[[seurat]], g2m.features=g2m.genes, s.features=s.genes)
      
      message("Finding Variable Genes...")
      
      seurat_filtered[[seurat]] <- Seurat::FindVariableFeatures(seurat_filtered[[seurat]], selection.method = "vst", nfeatures = n_variable_feats, verbose = F)
      
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  
  # integration # Order of the samples in the list may be important....so be careful
  
  seurat_anchors <- Seurat::FindIntegrationAnchors(object.list = seurat_filtered, dims = 1:npcs_integration)
  
  seurat_integrated <- Seurat::IntegrateData(anchorset = seurat_anchors, dims = 1:npcs_integration)
  
  # integrated analysis
  
  DefaultAssay(seurat_integrated) <- "integrated"
  
  set.seed(1993)
  
  message("Scaling Data Matrix...")
  seurat_integrated <- Seurat::ScaleData(seurat_integrated, verbose = FALSE)
  
  message("Performing PCA...")
  seurat_integrated <- Seurat::RunPCA(seurat_integrated, npcs = npcs_clustering, verbose = FALSE)
  
  message("Performing UMAP...")
  seurat_integrated <- Seurat::RunUMAP(seurat_integrated, reduction = "pca", dims = 1:npcs_clustering, seed.use=1993)
  
  message("Performing TSNE...")
  seurat_integrated <- Seurat::RunTSNE(seurat_integrated, reduction = "pca", dims = 1:npcs_clustering, seed.use=1993)
  
  message("Finding Neighbors...")
  seurat_integrated <- Seurat::FindNeighbors(seurat_integrated, reduction = "pca", dims = 1:npcs_clustering)
  
  message("Clustering...")
  seurat_integrated <- Seurat::FindClusters(seurat_integrated, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0))
  
  return(seurat_integrated)
}
