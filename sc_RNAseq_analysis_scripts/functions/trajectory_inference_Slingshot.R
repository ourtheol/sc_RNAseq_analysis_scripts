


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















## Identify temporally dynamic genes
# fit negative binomial GAM

sceGAM <- fitGAM(counts = counts(mut.traj),
                 sds=SlingshotDataSet(mut.traj))

## enter into the sceGAM object the RNA_snn_res.0.4 clustering
cl.0.4 <- colData(mut.traj)$RNA_snn_res.0.4
print(head(cl.0.4))
colData(sceGAM)$RNA_snn_res.0.4 <- cl.0.4

colnames(colData(sceGAM))     ## colData names(3): slingshot tradeSeq RNA_snn_res.0.4

# test for dynamic expression
ATres <- associationTest(sceGAM)

## access slingshot results: colData(sceGAM)$slingshot or sceGAM$slingshot
## access each curve:  print(colData(sceGAM)$slingshot$pseudotime.curve1)

topgenes <- rownames(ATres[order(ATres$pvalue), ])[1:100]
pst.ord <- order(as.numeric(as.character(colData(sceGAM)$slingshot$pseudotime.curve3)), decreasing=FALSE)    #######change curve
heatclus <- sceGAM$RNA_snn_res.0.4[pst.ord]


topgenes <- rownames(ATres[order(ATres$pvalue), ])[1:100]

# print(head(colData(sceGAM)$slingshot))  ### the dataframe of the curves

df <- colData(sceGAM)$slingshot
cell.ord <- rownames(df[order(as.numeric(as.character(df$pseudotime.curve3))),])  ###order the cells in increasing order according to pseudotime

# print(head(assays(sceGAM)$counts))
print(class(assays(sceGAM)$counts))
print(typeof(assays(sceGAM)$counts))
print(is.numeric(assays(sceGAM)$counts))

print(head(assays(sceGAM)$counts)[,1:4])

heatdata.cell.ord <- assays(sceGAM)$counts[, cell.ord]
print(head(heatdata.cell.ord)[,1:4])

heatdata.cell.ord.gene.ord <- heatdata.cell.ord[topgenes, ]
print(head(heatdata.cell.ord.gene.ord)[,1:4])

print(head(topgenes))
print(head(cell.ord))

heatclus <- sceGAM$RNA_snn_res.0.4[cell.ord]  ###edw ena ordering
print(head(heatclus))

print(head(sceGAM))
print(head(sceGAM$RNA_snn_res.0.4))

heatclus <- sceGAM$RNA_snn_res.0.4[cell.ord]
print(head(heatclus))


print(head(mut.traj$RNA_snn_res.0.4))
heatclus <- mut.traj$RNA_snn_res.0.4[cell.ord]
print(head(heatclus))


pst.ord <- order(as.numeric(as.character(colData(sceGAM)$slingshot$pseudotime.curve3)), decreasing=FALSE)
heatclus <- sceGAM$RNA_snn_res.0.4[pst.ord]
print(head(heatclus))


topgenes <- rownames(ATres[order(ATres$pvalue), ])[1:100]
pst.ord <- order(mut.traj$slingPseudotime_1, na.last = NA)    #######change curve
heatdata <- assays(sceGAM)$counts[topgenes, pst.ord]
print(head(heatdata))
heatclus <- mut.traj$RNA_snn_res.0.4[pst.ord]


purpleyellow <- c("#40004B", "#FFFF33")
pal <- colorRampPalette(purpleyellow)(100)


tiff(filename = paste0(out_path,"MYD88mut_traj_UMAP_heatmap_curve3_0.tiff"),width = 600, height = 600, units = "px")    #######change curve
print(heatmap(log1p(heatdata.cell.ord.gene.ord),
              col = pal,
              labCol = FALSE,
              ColSideColors = brewer.pal(11,'Spectral')[heatclus]))
# print(legend(x="topleft", legend="Expression", fill=pal))
dev.off()


tiff(filename = paste0(out_path,"MYD88mut_traj_UMAP_heatmap_curve3_1.tiff"),width = 600, height = 600, units = "px")    #######change curve
print(heatmap(log1p(heatdata.cell.ord.gene.ord),
              col = pal,
              labCol = FALSE,
              keep.dendro = TRUE,
              ColSideColors = brewer.pal(11,'Spectral')[heatclus]))
# print(legend(x="topleft", legend="Expression", fill=pal))
dev.off()


tiff(filename = paste0(out_path,"MYD88mut_traj_UMAP_heatmap_curve3_2.tiff"),width = 600, height = 600, units = "px")    #######change curve
print(heatmap(log1p(heatdata.cell.ord.gene.ord),
              col = pal,
              Rowv = NA, Colv = "Rowv",
              labCol = FALSE,
              ColSideColors = brewer.pal(11,'Spectral')[heatclus]))
# print(legend(x="topleft", legend="Expression", fill=pal))
dev.off()


tiff(filename = paste0(out_path,"MYD88mut_traj_UMAP_heatmap_curve3_3.tiff"),width = 600, height = 600, units = "px")    #######change curve
print(heatmap(log1p(heatdata.cell.ord.gene.ord),
              col = pal,
              Rowv = NA, Colv = "Rowv",
              labCol = TRUE,
              ColSideColors = brewer.pal(11,'Spectral')[heatclus]))
# print(legend(x="topleft", legend="Expression", fill=pal))
dev.off()


pdf(paste0(out_path,"MYD88mut_traj_UMAP_heatmap_curve3_3.pdf"))    #######change curve
print(heatmap(log1p(heatdata.cell.ord.gene.ord),
              col = pal,
              Rowv = NA, Colv = "Rowv",
              labCol = TRUE,
              ColSideColors = brewer.pal(11,'Spectral')[heatclus]))
# print(legend(x="topleft", legend="Expression", fill=pal))
dev.off()

png(filename = paste0(out_path,"MYD88mut_traj_UMAP_heatmap_curve3.png"),width = 600, height = 600, units = "px")    #######change curve
print(heatmap(log1p(heatdata),
              col = pal,
              Colv = NA,
              labCol = FALSE,
              ColSideColors = brewer.pal(11,'Spectral')[heatclus]))
# print(legend(x="topleft", legend="Expression", fill=pal))
dev.off()

#####################################




# Subset the sce object if needed
sce <- prep_seurat_for_slingshot(integrated_sobj)

mut.traj <- sce[, sce$condition == "MYD88mut"]

# Trajectory Analysis using Slingshot
mut.traj <- slingshot(mut.traj,
                      clusterLabels = "RNA_snn_res.0.4",
                      reducedDim = "UMAP", start.clus = 8)

P <- plot_slingshot_trajectories(mut.traj)
plot_slingshot_trajectories_per_lineage((mut.traj))




