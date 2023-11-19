library(Seurat)
library(SeuratData)
# library(sctransform)
# library(harmony)
library(ggplot2)
library(dplyr)
# library(grDevices)
library(cowplot)
library(patchwork)
library(purrr)
# library(tidyverse)

options(future.globals.maxSize = 4000 * 1024^6)

set.seed(1993)

select <- dplyr::select


# function: Read10x_CreateSeurat -----------------------------------------------
# Specify a path of the directory of the Cell Ranger filtered_feature_bc_matrix
# output and it will return a list all of the respective Seurat objects
# Note: Name the filtered_feature_bc_matrix directories accordingly to have nicely named Seurat objects that represent the biology of each sample

Read10x_CreateSeurat <- function(dir_path) { 
  
  set.seed(1993)
  
  require(Seurat)
  
  objects <- as.list(list.dirs(dir_path, full.names = FALSE, recursive=F))
  
  seurat_list <- list()
  
  for (object in objects) {
    
    print(object)
    
    seurat_data <- Read10X(data.dir = paste0(dir_path, object))
    
    seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                     min.features = 0,
                                     min.cells = 0,
                                     project = object)
    
    seurat_list[[object]] <- seurat_obj
    
    remove(seurat_data, seurat_obj)
  }
  
  return(seurat_list)
  
}


# function: MergeSeurat --------------------------------------------------------
# Merge all of the seurat objects in a list

MergeSeurat <- function(seurat_list) {

  require(Seurat)
  
  print(seurat_list)
  
  seurat_merged <- merge(x = seurat_list[[1]],
                         y = seurat_list[c(2:length(seurat_list))],
                         add.cell.id = names(seurat_list),
                         project = "Waldenstrom",
                         merge.data = F)#merge.data=TRUE will keep normalized data matrices and raw
  
  seurat_merged[["percent.mt"]] <- PercentageFeatureSet(seurat_merged, pattern = "^MT-") # calculate mitochondrial gene percentage per cell
  
  seurat_merged[["percent.rb"]] <- PercentageFeatureSet(seurat_merged, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA|^RACK1") # calculate ribosomal gene percentage per cell
  
  return(seurat_merged)
  
}


# function: FilterSeurat -------------------------------------------------------
# Filter all cells in a merged Seurat object with universal filtering parameters 

FilterSeurat <- function(seurat_merged, 
                         min_genes = 200, 
                         max_genes = 2500, 
                         max_mito_genes = 15, 
                         min_non_zero_cells_per_gene = 10) {
						 
  require(Seurat)
  require(dplyr)
  
  # filter cells
  message("cell and gene numbers before filtering")
  print(dim(seurat_merged))
  
  seurat_filtered <- subset(x = seurat_merged, subset = nFeature_RNA > min_genes & nFeature_RNA < max_genes & percent.mt < max_mito_genes)
  
  message("cell and gene numbers after cell filtering")
  print(dim(seurat_filtered))
  
  # filter genes
  
  counts <- GetAssayData(object = seurat_filtered, slot = "counts")# Output a logical vector for every gene on whether the more than zero counts per cell
  
  nonzero <- (counts > 0)
  
  keep_genes <- Matrix::rowSums(nonzero) >= min_non_zero_cells_per_gene # Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
  
  filtered_counts <- counts[keep_genes, ] # Only keeping those genes expressed in more than 10 cells
  
  seurat_filtered <- CreateSeuratObject(filtered_counts, meta.data = seurat_filtered@meta.data)# Reassign to filtered seurat object
  
  message("cell and gene numbers after gene filtering")
  print(dim(seurat_filtered))
  
  return(seurat_filtered)
}


# function: PlotQC -------------------------------------------------------------
# Create .pdf plots of metrics used for filtering before and after

PlotQC <- function(seurat_merged, 
                   seurat_filtered, 
                   out_path = out_path) {
  
  # before filtering
  
  pdf(paste0(out_path,"qc_before_filtering.pdf"), width = 150, height = 40, bg = "white")
  print(VlnPlot(seurat_merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), group.by = "orig.ident", ncol = 4)&theme(aspect.ratio = 1))
  dev.off()
  
  # after filtering
  
  pdf(paste0(data_path,"qc_after_filtering.pdf"), width = 150, height = 40, bg = "white")
  print(VlnPlot(seurat_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), group.by = "orig.ident", ncol = 4)&theme(aspect.ratio = 1))
  dev.off()
   
}


# function: SeparateSeuratAnalysis ---------------------------------------------
# Standard 'Seurat' analysis and clustering of cells in each Seurat object or merged Seurat objects

SeparateSeuratAnalysis <- function(seurat_filtered,
                                   n_variable_feats, #n_variable_feats between 2000-3000 recommended
                                   clustering_algorithm,
                                   clustering_resolution,
                                   out_path) {

  require(Seurat)
  require(scran)
  require(dplyr)
  
  # A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  
  # We can segregate this list into markers of G2/M phase and markers of S phase
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes

  if (!is.list(seurat_filtered)) {

    seurat_filtered <- Seurat::SplitObject(seurat_filtered, split.by = "orig.ident")

  }

  select <- dplyr::select

  #load("C:/Users/Fotini/Documents/Meduoa_research/scrna_seq/data/cycle.rda")#load cell cycle genes must be downloaded and loaded

  for (seurat in 1:length(seurat_filtered)) {

    tryCatch({

      set.seed(1993)

      print("Analysis of sample: ")

      print(seurat_filtered[[seurat]]$orig.ident[[1]])

      print("Normalizing...")

      seurat_filtered[[seurat]] <- Seurat::NormalizeData(seurat_filtered[[seurat]])#, assay = NULL, normalization.method = "LogNormalize", scale.factor = 10000, margin = 1, verbose = F)

      print("Cell Cycle Scoring...")

      seurat_filtered[[seurat]] <- Seurat::CellCycleScoring(seurat_filtered[[seurat]], g2m.features=g2m.genes, s.features=s.genes)

      print("Finding Variable Genes...")

      seurat_filtered[[seurat]] <- Seurat::FindVariableFeatures(seurat_filtered[[seurat]], selection.method = "vst", nfeatures = n_variable_feats, verbose = F)

      print("Scaling Data Matrix...")

      seurat_filtered[[seurat]] <- Seurat::ScaleData(seurat_filtered[[seurat]], features = rownames(seurat_filtered[[seurat]]), verbose = F)#added vars.to.regress, vars.to.regress = c("S.Score", "G2M.Score")

      
      print("Performing PCA...")
      
      seurat_filtered[[seurat]] <- Seurat::RunPCA(object = seurat_filtered[[seurat]], features = VariableFeatures(object = seurat_filtered[[seurat]]), do.print = F, npcs = 40, verbose = F, seed.use = 1993)
      
      # Determine how many PCs capture the majority of the variation in the data
      # Elbow plot
      pdf(paste0(out_path,seurat_filtered[[seurat]]$orig.ident[[1]],"elbowplot.pdf"), width = 15, height = 10, bg = "white")
      print(ElbowPlot(seurat_filtered[[seurat]],ndims = 40))
      dev.off()
      # Determine percent of variation associated with each PC
      pct <- seurat_filtered[[seurat]][["pca"]]@stdev / sum(seurat_filtered[[seurat]][["pca"]]@stdev) * 100
      # Calculate cumulative percents for each PC
      cumu <- cumsum(pct)
      # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
      co1 <- which(cumu > 90 & pct < 5)[1]
      # Determine the difference between variation of PC and subsequent PC
      co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
      # last point where change of % of variation is more than 0.1%.
      # Minimum of the two calculation
      pcs <- min(co1, co2)
      message(paste0(seurat_filtered[[seurat]]$orig.ident[[1]], "number of PCs ", pcs))
            
      message("Executing UMAP Algorithm...")

      seurat_filtered[[seurat]] <- Seurat::RunUMAP(seurat_filtered[[seurat]], reduction = "pca", dims = 1:pcs, verbose = F, seed.use = 1993)

      message("Executing TSNE Algorithm...")

      seurat_filtered[[seurat]] <- Seurat::RunTSNE(seurat_filtered[[seurat]], reduction = "pca", dims = 1:pcs, verbose = F, seed.use = 1993)

      message("Finding Neighbors...")

      seurat_filtered[[seurat]] <- Seurat::FindNeighbors(seurat_filtered[[seurat]], dims = 1:pcs, verbose = F)

      message("Clustering...")

      seurat_filtered[[seurat]] <- Seurat::FindClusters(seurat_filtered[[seurat]], algorithm = clustering_algorithm, resolution = clustering_resolution, verbose = F)


    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }

  return(seurat_filtered)
}


# function: DifferentiallyExpressedGenes ---------------------------------------
# Standard 'Seurat' FindAllMarkers wrapper function for the identification of
# differentially expressed genes between clusters using the clustering result of preference (the appropriate Ident)
# This can be applied in a single Seurat object or a list of Seurat objects.
# The output of this function is a list of lists with the results for each Seurat object.

DifferentiallyExpressedGenes <- function(seurat_analyzed, clustering) { #list of Seurat objects or single Seurat object #the clustering ident to use for the de analysis

  require(Seurat)
  require(scran)
  require(dplyr)

  select <- dplyr::select

  if (!is.list(seurat_analyzed)) {

    message("Converting into a list...")

    temp <- list()

    temp[[print(as.character(substitute(seurat_analyzed)))]] <- seurat_analyzed

    seurat_analyzed <- temp

  }

  de_analysis_results <- list()

  for (seurat in 1:length(seurat_analyzed)) {

    tryCatch({

      message("Creating a list to store all of the results:")

      de_analysis_results[[print(names(seurat_analyzed[seurat]))]][["seurat_object"]] <- seurat_analyzed[[seurat]]

      # count cells in each cluster of cells

      message("Counting cells in each cluster of cells:")

      Idents(object = de_analysis_results[[seurat]]$seurat_object) <- clustering

      de_analysis_results[[seurat]][["clustering_cell_counts"]] <- FetchData(de_analysis_results[[seurat]]$seurat_object, vars = c("ident", "orig.ident")) %>%
        dplyr::count(ident, orig.ident) %>%
        tidyr::spread(ident, n)

      # identify marker genes in each cluster

      message("Identifying marker genes in each cluster:")

      DefaultAssay(de_analysis_results[[seurat]]$seurat_object) <- "RNA"

      de_analysis_results[[seurat]][["clustering_markers"]] <- FindAllMarkers(de_analysis_results[[seurat]]$seurat_object,
                                                                              assay = "RNA",
                                                                              only.pos = T,
                                                                              min.pct = 0.2,
                                                                              logfc.threshold = 0.5,
                                                                              verbose = F)
      # top 10 markers of each cluster

      message("Selecting the top 10 marker genes in each cluster:")

      de_analysis_results[[seurat]][["clustering_markers_top10"]] <- de_analysis_results[[seurat]][["clustering_markers"]] %>%
        group_by(cluster) %>%
        top_n(n = 10, wt = avg_log2FC)#choosing top highest expressed markers with the highest log2FC in each cluster

      de_analysis_results[[seurat]][["clustering_markers_top10"]] <- unique(de_analysis_results[[seurat]][["clustering_markers_top10"]]$gene)


    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }

  return(de_analysis_results)

}


# function: PrintPlots ---------------------------------------------------------
## Plotting of some basic plots for the analysis

PrintPlots <- function(de_analysis_results, out_path) {

  require(Seurat)
  require(scran)
  require(ggplot2)
  require(ggrepel)
  require(ggnewscale)
  require(ggbeeswarm)
  require(ggthemes)
  require(scales)
  require(RColorBrewer)
  require(Polychrome)
  require(viridis)
  require(ggpubr)
  require(reshape)
  require(cowplot)
  require(patchwork)
  require(purrr)
  # require(tidyverse)

  for (seurat in 1:length(de_analysis_results)) {

    tryCatch({

      # create all the plots first

      message("Making the plots...")
      
      message(paste0("Plots for " ,de_analysis_results[[seurat]]$orig.ident[[1]]))

      message("Creating UMAP...")
      umap <- DimPlot(de_analysis_results[[seurat]]$seurat_object, reduction = "umap", label = T, label.size = 4)&theme(aspect.ratio=1)
      Idents(de_analysis_results[[seurat]]$seurat_object) <- "Phase"
      umap.2 <- DimPlot(de_analysis_results[[seurat]]$seurat_object, reduction = "umap", label = F, label.size = 4)&theme(aspect.ratio=1)
      Idents(de_analysis_results[[seurat]]$seurat_object) <- "seurat_clusters"
      
      message("Creating TSNE...")
      tsne <- DimPlot(de_analysis_results[[seurat]]$seurat_object, reduction = "tsne", label = T, label.size = 4)&theme(aspect.ratio=1)
      # message("Creating PCA...")
      # pca <- DimPlot(de_analysis_results[[seurat]]$seurat_object, reduction = "pca", label = T, label.size = 4)&theme(aspect.ratio=1)
      
      message("Plot variable features...")
      # Identify the 15 most highly variable genes
      top15 <- head(VariableFeatures(de_analysis_results[[seurat]]$seurat_object), 15)
      # plot variable features with and without labels
      plot1 <- VariableFeaturePlot(de_analysis_results[[seurat]]$seurat_object)
      plot2 <- LabelPoints(plot = plot1, points = top15, repel = TRUE)
      
      message("Creating DOTPLOT...")
      dotplot <- DotPlot(de_analysis_results[[seurat]]$seurat_object,
                         assay = "RNA",
                         features = de_analysis_results[[seurat]][["clustering_markers_top10"]],
                         cluster.idents = F,
                         cols = c("white", "red"),
                         scale.by = "radius") + RotatedAxis() + theme(aspect.ratio = 0.15,
                                                                      legend.text=element_text(size=8),
                                                                      legend.title=element_text(size=8))
      message("Creating HEATMAP...")
      heatmap <- DoHeatmap(de_analysis_results[[seurat]]$seurat_object,
			   assay = "RNA",
                           features = de_analysis_results[[seurat]][["clustering_markers_top10"]])+ theme(aspect.ratio = 2,
                                                                                                          legend.text=element_text(size=8),
                                                                                                          legend.title=element_text(size=8))
      message("Creating FEATUREPLOTS...")
      featureplots <- FeaturePlot(de_analysis_results[[seurat]]$seurat_object, reduction = "umap", features = de_analysis_results[[seurat]][["clustering_markers_top10"]], label = T, label.size = 2, ncol = 10)&theme(aspect.ratio=1)

      message("Printing plots....")

      # print plots in PDF
      pdf(paste0(out_path, de_analysis_results[[seurat]]$seurat_object$orig.ident[[1]], "_results_plots.pdf"), width = 30, height = 20)

      # umap
      print(umap)
      
      print(umap.2)

      # tsne
      print(tsne)
      
      # variable features
      print(plot1 + plot2)

      # # pca
      # print(pca)

      # dotplot of top 10 markers in each cluster
      print(dotplot)

      # heatmap
      print(heatmap)

      # featureplots of top markers
      print(featureplots)

      dev.off()

      write.table(de_analysis_results[[seurat]][["clustering_markers"]], file = paste0(out_path,  de_analysis_results[[seurat]]$seurat_object$orig.ident[1], "_clustering_markers.txt"), row.names=FALSE, sep="\t", quote = FALSE )

      write.table(de_analysis_results[[seurat]][["clustering_markers_top10"]], file = paste0(out_path,  de_analysis_results[[seurat]]$seurat_object$orig.ident[1], "_clustering_markers_top10.txt"), row.names=FALSE,sep="\t", quote = FALSE )

      remove(umap, tsne, pca, dotplot, heatmap, featureplots) 

    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }

}
