library("Seurat")
library("slingshot")
library("RColorBrewer")
library("tradeSeq")


# function: prep_seurat_for_slingshot ------------------------------------------
# insert a seurat object in order to be turned into an sce object for input in slingshot

prep_seurat_for_slingshot <- function(Sobj){
  
  # Working with sce objects: insert seurat object and turn it into sce object
  sce <- as.SingleCellExperiment(Sobj)
  
  # Gene filtering
  # filter genes down to potential cell-type markers
  # at least M (15) reads in at least N (15) cells
  geneFilter <- apply(assays(sce)$counts,1,function(x){
    sum(x >= 3) >= 10
  })
  sce <- sce[geneFilter, ]
  
  return(sce)
}


# function: plot_slingshot_trajectories ----------------------------------------
# insert a slingshot object and get figures of the trajectories

plot_slingshot_trajectories <- function(slingshot_obj){
  
  colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
  plotcol <- colors[cut(slingshot_obj$slingPseudotime_1, breaks=100)]
  
  p <- plot(reducedDims(slingshot_obj)$UMAP, col = plotcol, pch=16, asp = 1)
  P <- p + lines(SlingshotDataSet(slingshot_obj), lwd=2, col='black')
  
  tiff(filename = "trajectories_clusters.tiff")
  plot(as.data.frame(reducedDims(slingshot_obj)$UMAP), col = plotcol, pch=16, asp = 1, xlab="UMAP_1", ylab="UMAP_2") +
    lines(SlingshotDataSet(slingshot_obj), lwd=2, type = 'lineages', col='black')&theme(aspect.ratio=1)
  dev.off()

  tiff(filename = "trajectories.tiff")
  plot(as.data.frame(reducedDims(slingshot_obj)$UMAP), col = plotcol, pch=16, asp = 1, xlab="UMAP_1", ylab="UMAP_2")+
    lines(SlingshotDataSet(slingshot_obj), lwd=2, col='black')&theme(aspect.ratio=1)
  dev.off()
  
  return(P)
}


# function: plot_slingshot_trajectories_per_lineage ----------------------------
# Which cells are in which lineage: plot the pseudotime values for each lineage.

plot_slingshot_trajectories_per_lineage <- function(slingshot_obj){
  
  nc <- 3
  pt <- slingPseudotime(slingshot_obj)
  nms <- colnames(pt)
  nr <- ceiling(length(nms)/nc)
  pal <- viridis(100, end = 0.95)
  par(mfrow = c(nr, nc))
  
  for (i in nms) {
    colors <- colorRampPalette(brewer.pal(11,'Spectral')[-10])(100)
    plotcol <- colors[cut(pt[,i], breaks=100)]
    tiff(filename =paste0("pseudotimeforeachlineage_clusters_lineage_", i, ".tiff"))
    plot(as.data.frame(reducedDims(slingshot_obj)$UMAP), col = plotcol, pch=16, main = i, asp = 1, xlab="UMAP_1", ylab="UMAP_2") +
      lines(SlingshotDataSet(slingshot_obj), lwd=2, type = 'lineages', col='black')&theme(aspect.ratio=1)
    dev.off()
    
    tiff(filename =paste0("pseudotimeforeachlineage_lineage_", i, ".tiff"))
    plot(as.data.frame(reducedDims(slingshot_obj)$UMAP), col = plotcol, pch=16, main = i, asp = 1, xlab="UMAP_1", ylab="UMAP_2") +
      lines(SlingshotDataSet(slingshot_obj), lwd=2, col='black')&theme(aspect.ratio=1)
    dev.off()
  }
}


# function: plot_slingshot_heatmap ---------------------------------------------
# Identify temporally dynamic genes and plot heatmap
# insert a slingshot object, the clustering used in this analysis, the selected curve genes to be plotted (eg pseudotime.curve1, pseudotime.curve2 etc)

plot_slingshot_trajectories_per_lineage <- function(slingshot_obj, clustering, curve){
  
  # fit negative binomial GAM
  
  sceGAM <- fitGAM(counts = counts(slingshot_obj),
                   sds=SlingshotDataSet(slingshot_obj))
  
  ## enter into the sceGAM object the clustering column name
  clust <- colData(slingshot_obj)$clustering
  colData(sceGAM)$clustering <- clust
  
  # test for dynamic expression
  ATres <- associationTest(sceGAM)

  topgenes <- rownames(ATres[order(ATres$pvalue), ])[1:100]
  pst.ord <- order(colData(sceGAM)$slingshot$curve, na.last = NA)    
  heatdata <- assays(sceGAM)$counts[topgenes, pst.ord]
  heatclus <- sceGAM$clustering[pst.ord]
  
  tiff(filename = paste0("heatmap_", curve,".tiff"), width = 600, height = 600, units = "px")    
  print(heatmap(log1p(heatdata),
                Colv = NA,
                labCol = FALSE,
                ColSideColors = brewer.pal(11,'Spectral')[heatclus]))
  dev.off()
}



# # Example of using the above functions:
# 
# # Subset the sce object if needed
# sce <- prep_seurat_for_slingshot(integrated_sobj)
# 
# mut.traj <- sce[, sce$condition == "MYD88mut"]
# 
# # Trajectory Analysis using Slingshot
# mut.traj <- slingshot(mut.traj,
#                       clusterLabels = "RNA_snn_res.0.4",
#                       reducedDim = "UMAP", start.clus = 8)
# 
# P <- plot_slingshot_trajectories(mut.traj)
# plot_slingshot_trajectories_per_lineage(mut.traj)
# 
# plot_slingshot_trajectories_per_lineage(mut.traj, seurat_clusters, pseudotime.curve1)


