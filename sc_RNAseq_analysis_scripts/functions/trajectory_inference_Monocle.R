library("Seurat")
library("ggplot2")
library("tidyverse")
library("magrittr")
library("monocle")
library("pheatmap")
library("RColorBrewer")
library("colorRamps")


# function: prep_seurat_for_monocle ------------------------------------------
# insert a seurat object and select the clustering you want to infer a trajectory for
# returns a ready to be analyzed with monocle object

prep_seurat_for_slingshot <- function(data_integrated.seurat_object, selected_clustering){
  
  # Extract data, phenotype data, and feature data from the SeuratObject manually and build the monocle CDS
  
  # Extract the read count matrix from Seurat object, and then create a new Monocle object
  # sparse matrix: cols --> cells, rows --> genes
  data <- as(as.matrix(GetAssayData(object = data_integrated.seurat_object, slot = "counts")), 'sparseMatrix')
  
  pd <- new('AnnotatedDataFrame', data = data_integrated.seurat_object@meta.data)
  
  # Rename the selected clustering column as "Cluster" so that monocle thinks clustering has been performed
  colnames(pd)[colnames(pd)==selected_clustering] <- "Cluster" 
  
  fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
  
  fd <- new('AnnotatedDataFrame', data = fData)
  
  
  #Construct monocle cds
  # Choosing a distribution for your data
  # UMIs or read counts are better modeled with the negative binomial
  data_integratedmonocle <- newCellDataSet(data,
                                           phenoData = pd,
                                           featureData = fd,
                                           lowerDetectionLimit = 0.5,
                                           expressionFamily = negbinomial.size())
  
  
  # Extract the UMAP matrix calculated by Seurat, and then import it to monocle pipeline
  umap.embedding <- Embeddings(object = data_integrated.seurat_object, reduction = "umap")
  colnames(umap.embedding)<-NULL   ## remove the column names UMAP_1 and UMAP_2
  data_integratedmonocle@reducedDimA <- t(umap.embedding)
  
  ## All the above because plot_cell_clusters() checks the following:
  ## if (is.null(cds@reducedDimA) | length(pData(cds)$Cluster) == 0){
  ##   stop("Error: Clustering is not performed yet. Please call clusterCells() before calling this function.")
  # print(head(data_integratedmonocle@reducedDimA))
  # print(length(pData(data_integratedmonocle)$Cluster))
  
  # plotting to check that the Seurat clusters have been correctly imported
  p <- plot_cell_clusters(data_integratedmonocle, x = 1, y = 2, color_by = "Cluster", verbose = T)&theme(aspect.ratio=1)
  pdf("monocle_integrated_clustering_onUMAP.pdf", width = 10, height = 5, bg = "white")
  print(p)
  dev.off()
  
  return(data_integratedmonocle)
}


# function: monocle_trajectories  ----------------------------------------------
# insert a ready to be analyzed with monocle object and get figures of the trajectories

monocle_trajectories <- function(HSMM){
  
  HSMM <- estimateSizeFactors(data_integratedmonocle)
  HSMM <- estimateDispersions(HSMM)
  
  ## Filtering low-quality cells --> num_cells_expressed col into the fData df
  HSMM <- detectGenes(HSMM, min_expr = 0.1)
  
  ## genes expressed in > 50 cells
  HSMM_expressed_genes <-  row.names(subset(fData(HSMM), num_cells_expressed >= 50))
  
  ## After we confirm the clustering makes sense,
  ## perform differential gene expression test as a way to extract the genes that distinguish them.
  clustering_DEG_genes <- differentialGeneTest(HSMM[HSMM_expressed_genes,],
                                               fullModelFormulaStr = '~Cluster',
                                               cores = 1)
  
  ## select the top 1000 significant genes as the ordering genes.
  HSMM_ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]
  
  # Once we have a list of gene ids to be used for ordering, we need to set them in the HSMM object,
  # because the next several functions will depend on them.
  HSMM <- setOrderingFilter(HSMM, ordering_genes = HSMM_ordering_genes)
  pdf("plot_ordering_genes.pdf", width = 10, height = 5, bg = "white")
  print(plot_ordering_genes(HSMM))
  dev.off()
  
  # reduce the space down to one with two dimensions
  HSMM <- reduceDimension(HSMM, max_components = 2, method = 'DDRTree')
  
  # Now that the space is reduced, order cells along the trajectory,
  # we often have to call orderCells again using the root_state argument to specify the beginning, check next function
  HSMM <- orderCells(HSMM)
  # HSMM <- orderCells(HSMM, root_state = 6)

  # plot_complex_cell_trajectory() gets the x and y axis from @reducedDimS
  pdf(paste0(out_path,"monocle_integrated_traj_colorbyClusterState.pdf"), width = 10, height = 10, bg = "white")
  print(plot_cell_trajectory(HSMM, color_by = "State")&theme(aspect.ratio=1))
  print(plot_cell_trajectory(HSMM, color_by = 'as.factor(Cluster)')&theme(aspect.ratio=1)) ##, color_by = 'as.factor(Cluster)'
  print(plot_complex_cell_trajectory(HSMM, color_by = "State")&theme(aspect.ratio=1))
  dev.off()
  
  ## "facet" the trajectory plot so it's easier to see where each of the states are located:
  pdf(paste0(out_path,"monocle_integrated_traj_facet.pdf"), width = 100, height = 10, bg = "white")
  print(plot_cell_trajectory(HSMM, color_by = "State") + facet_wrap(~State, nrow = 1)&theme(aspect.ratio=1))
  print(plot_cell_trajectory(HSMM, color_by = 'as.factor(Cluster)') + facet_wrap(~Cluster, nrow = 1)&theme(aspect.ratio=1))
  dev.off()
  
  return(HSMM)
}


# function: GM_state  ----------------------------------------------------------
# "State" --> Monocle's term for the segment of the tree
# The function below is handy for identifying the State which contains most of the cells from a selected cluster.
# We can then pass that to orderCells

GM_state <- function(cds, metadatacol, colID){
  if (length(unique(pData(cds)$metadatacol)) > 1){
    T0_counts <- table(pData(cds)$metadatacol, pData(cds)$metadatacol)[,colID]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}



# # Example of using the above functions:
# 
# data_integratedmonocle <- prep_seurat_for_slingshot(data_integrated.seurat_object, "seurat_clusters")
# HSMM <- monocle_trajectories(data_integratedmonocle)
# HSMM <- orderCells(HSMM, root_state = GM_state(HSMM, Cluster, "8"))
# 
# pdf(paste0(out_path,"monocle_integrated_traj_colorbyPseudotime_fromClus8.pdf"), width = 10, height = 5, bg = "white")
# print(plot_cell_trajectory(HSMM, color_by = "Pseudotime"))
# dev.off()
# 
# 
# # Test each gene for differential expression as a function of pseudotime or according to other covariates as specified.
# diff_test_res <- differentialGeneTest(HSMM,
#                                       fullModelFormulaStr = "~sm.ns(Pseudotime)",
#                                       cores = 8)
# 
# # sort the diff_test_res according to pval (3rd col) in increasing order
# # and keep the 100 most statistically significant genes
# diff_test_res_rearranged <- diff_test_res[order(diff_test_res[,3],decreasing=FALSE),][1:100,]
# # make sure that the genes plotted have qval < 0.1
# sig_gene_names <- row.names(subset(diff_test_res_rearranged, qval < 0.1))
# 
# 
# 
# # Heatmap plotting
# 
# # plot_pseudotime_heatmap() plots ONLY 100 cells on the y axis
# # use newdata to print the pseudotime values as in plot_pseudotime_heatmap()
# newdata <- data.frame(Pseudotime = seq(min(pData(HSMM[sig_gene_names,])$Pseudotime), max(pData(HSMM[sig_gene_names,])$Pseudotime),length.out = 100))
# 
# # get the states
# binner <- function(cds_object){
#   df <- data.frame(pData(cds_object))
#   df <- df[,c("Pseudotime", "State")]
#   df <- df[order(df$Pseudotime, decreasing = F),]
#   len <- length(df$Pseudotime)
#   bin<-round(len/100)
#   State <- c()
#   value <- c()
#   for(i in 0:99){
#     if(i < 99){
#       start <- 1+(bin*i)
#       stop <- bin+(bin*i)
#       value <- median(as.numeric(as.vector(df$State[c(start:stop)])))
#       State <- c(State, value)
#     }
#     else{
#       State <- c(State, value)
#     }
#   }
#   return(as.data.frame(State))
# }
# 
# bin <- binner(HSMM[sig_gene_names,])
# 
# # transform the numerical values of States into categorical values for the legend
# bin$State <- as.factor(bin$State)
# 
# annotation_col = data.frame(
#   # State = bin,
#   Pseudotime = newdata
# )
# 
# # use the color palette Blues (white -> dark blue) and reverse the colors (dark blue -> light blue)
# blues.pseudotime <- rev(colorRampPalette(brewer.pal(9, "Blues"))(100))[1:60]
# ann_colors = list(
#   Pseudotime = blues.pseudotime)  #,
#   # State = c("1"="#F8766D", "2"="#CD9600", "3"="#A3A500", "4"="#53B400", "5"="#00BE67", "6"="#00BFC4", "7"="#00A9FF", "8"="#A58AFF", "9"="#FB61D7")
# # )
# 
# 
# 
# source("functions/plot_pseudotime_heatmap_monocle_adjustments.R")   # check for bugs
# 
# heatmap <- plot_pseudotime_heatmap_2(HSMM[sig_gene_names,],
#                                    cluster_rows = T,
#                                    num_clusters = 5,
#                                    add_annotation_col = annotation_col,
#                                    annotation_colors = ann_colors,
#                                    cores = 8,
#                                    show_rownames = T,
#                                    return_heatmap = T)
