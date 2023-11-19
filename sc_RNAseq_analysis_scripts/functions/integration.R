library("Seurat")
library("harmony")



# function: harmony.integration ------------------------------------------------
# insert a list of seurat objects, 
# seurat objects will be merged and integrated using harmony

harmony.integration <- function(Sobj_list,
                                n_variable_feats, #n_variable_feats between 2000-3000 recommended
                                clustering_algorithm,
                                clustering_resolution) {
  
  require("Seurat")
  require("harmony")
  
  # Merge seurat objects
  seurat_merged <- merge(x = Sobj_list[[1]],
                         y = Sobj_list[c(2:length(Sobj_list))],
                         add.cell.id = names(Sobj_list),
                         project = "Waldenstrom",
                         merge.data = F)#merge.data=TRUE will keep normalized data matrices and raw
  
  
  harmony_integrated <- NormalizeData(object = seurat_merged)
  
  
  # A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  
  # We can segregate this list into markers of G2/M phase and markers of S phase
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  
  
  harmony_integrated <- CellCycleScoring(harmony_integrated, g2m.features=g2m.genes, s.features=s.genes)
  
  
  harmony_integrated <- FindVariableFeatures(object = harmony_integrated, selection.method = "vst", nfeatures = n_variable_feats, verbose = T)
  
  
  harmony_integrated <- ScaleData(harmony_integrated)#, vars.to.regress = c("S.Score", "G2M.Score"))
  
  
  # Determine how many PCs capture the majority of the variation in the data
  # Elbow plot
  pdf("harmony_integrated_elbowplot.pdf", width = 15, height = 10, bg = "white")
  print(ElbowPlot(harmony_integrated,ndims = 40))
  dev.off()
  # Determine percent of variation associated with each PC
  pct <- harmony_integrated[["pca"]]@stdev / sum(harmony_integrated[["pca"]]@stdev) * 100
  # Calculate cumulative percents for each PC
  cumu <- cumsum(pct)
  # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
  co1 <- which(cumu > 90 & pct < 5)[1]
  # Determine the difference between variation of PC and subsequent PC
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  # last point where change of % of variation is more than 0.1%.
  # Minimum of the two calculation
  pcs <- min(co1, co2)
  message(paste0("number of PCs ", pcs))
  
  harmony_integrated <- RunPCA(object = harmony_integrated, do.print = F, npcs = pcs, verbose = T, seed.use = 1993)
  
  
  harmony_integrated <- RunHarmony(object = harmony_integrated, assay.use = "RNA", group.by.vars = "orig.ident", reduction = "pca", dims.use = 1:pcs, verbose = T, seed.use = 1993)
  
  message("Executing UMAP Algorithm...")
  harmony_integrated <- RunUMAP(harmony_integrated, reduction = "harmony", dims = 1:pcs, seed.use = 1993)
  
  message("Finding Neighbors...")
  harmony_integrated <- FindNeighbors(harmony_integrated, reduction = "harmony", dims = 1:pcs) 
    
  message("Clustering...")
  harmony_integrated <- Seurat::FindClusters(harmony_integrated, algorithm = clustering_algorithm, resolution = clustering_resolution, verbose = F)
  
  
  return(harmony_integrated) 
}




# function: seurat.integration -------------------------------------------------
# Standard 'Seurat' integration and analysis of cells in a list of Seurat objects

seurat.integration <- function(Sobj_list,
                               npcs_integration,
                               clustering_algorithm,
                               clustering_resolution) {
  
  require("Seurat")
  require("scran")
  require("dplyr")
  
  # integration # Order of the samples in the list may be important....so be careful
  
  seurat_anchors <- Seurat::FindIntegrationAnchors(object.list = Sobj_list, dims = 1:npcs_integration)
  
  seurat_integrated <- Seurat::IntegrateData(anchorset = seurat_anchors, dims = 1:npcs_integration)
  
  # integrated analysis
  
  DefaultAssay(seurat_integrated) <- "integrated"
  
  set.seed(1993)
  
  message("Scaling Data Matrix...")
  seurat_integrated <- Seurat::ScaleData(seurat_integrated, verbose = FALSE)
  
  
  # Determine how many PCs capture the majority of the variation in the data
  # Elbow plot
  pdf("harmony_integrated_elbowplot.pdf", width = 15, height = 10, bg = "white")
  print(ElbowPlot(harmony_integrated,ndims = 40))
  dev.off()
  # Determine percent of variation associated with each PC
  pct <- harmony_integrated[["pca"]]@stdev / sum(harmony_integrated[["pca"]]@stdev) * 100
  # Calculate cumulative percents for each PC
  cumu <- cumsum(pct)
  # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
  co1 <- which(cumu > 90 & pct < 5)[1]
  # Determine the difference between variation of PC and subsequent PC
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  # last point where change of % of variation is more than 0.1%.
  # Minimum of the two calculation
  pcs <- min(co1, co2)
  message(paste0("number of PCs ", pcs))
  
  message("Performing PCA...")
  seurat_integrated <- Seurat::RunPCA(seurat_integrated, npcs = pcs, verbose = FALSE)
  
  message("Performing UMAP...")
  seurat_integrated <- Seurat::RunUMAP(seurat_integrated, reduction = "pca", dims = 1:pcs, seed.use=1993)
  
  message("Performing TSNE...")
  seurat_integrated <- Seurat::RunTSNE(seurat_integrated, reduction = "pca", dims = 1:pcs, seed.use=1993)
  
  message("Finding Neighbors...")
  seurat_integrated <- Seurat::FindNeighbors(seurat_integrated, reduction = "pca", dims = 1:pcs)
  
  message("Clustering...")
  seurat_integrated <- Seurat::FindClusters(seurat_integrated, algorithm = clustering_algorithm, resolution = clustering_resolution, verbose = F)
  
  return(seurat_integrated)
}
