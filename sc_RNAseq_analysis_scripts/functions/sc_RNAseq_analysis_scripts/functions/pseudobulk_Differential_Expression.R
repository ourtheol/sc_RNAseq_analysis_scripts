library("Seurat")
library("ggplot2")
library("ggrepel")


# function: pseudobulk_differential_expession ----------------------------------
# insert a seurat object, select the column name according to each pseudobulking is performed (condition),
# genes to be ignored, identity class 1, identity class 2 for comparison
# returns the differentially expressed genes

pseudobulk_differential_expession <- function(Sobj, condition, remove_genes, ident_1, ident_2){
  
  # Set identity classes to an existing column in meta.data
  Idents(Sobj) <- condition
  
  # Averaged feature expression by identity class
  avg.patient.cells <- log1p(AverageExpression(Sobj, verbose = FALSE)$RNA)
  avg.patient.cells.df <- data.frame(avg.patient.cells)
  
  DE.genes <- FindMarkers(Sobj,
                          features = rownames(Sobj@assays$RNA)[! rownames(Sobj@assays$RNA) %in% remove_genes],
                          ident.1 = ident_1, 
                          ident.2 = ident_2, 
                          verbose = FALSE) # there is default value: logfc.threshold = 0.25

  
  # The significantly differentially expressed genes:
  # Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange positive or negative respectively)
  DE.genes$diffexpressed <- "NO"
  # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"
  DE.genes$diffexpressed[DE.genes$avg_log2FC > 0.6 & DE.genes$p_val < 0.05] <- "UP"
  # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
  DE.genes$diffexpressed[DE.genes$avg_log2FC < -0.6 & DE.genes$p_val < 0.05] <- "DOWN"
  
  # Create a new column "delabel" to de, that will contain the name of the differentially expressed genes (NA in case they are not)
  DE.genes$delabel <- NA
  DE.genes$delabel[DE.genes$diffexpressed != "NO"] <- rownames(DE.genes)[DE.genes$diffexpressed != "NO"]
  
  ## write in a table the DE (UP, DOWN) genes
  DE.genes.to.print <- subset(DE.genes, diffexpressed != "NO")
  write.csv(x = DE.genes.to.print, file = "differentially_expressed_genes.csv")

  return(DE.genes.to.print)
}


# function: scatter_plot_DEgenes -----------------------------------------------
# insert a seurat object, the column name according to each pseudobulking is performed (condition),
# which conditions to be plotted on the x, y axis, the differentially expressed genes (genes_to_label),
# how many of the differentially expressed genes to be labeled on the plot (top_genes)

scatter_plot_DEgenes <- function(Sobj, condition, x_axis, y_axis, genes_to_label, top_genes){
  
  # Set identity classes to an existing column in meta.data
  Idents(Sobj) <- condition
  
  # Averaged feature expression by identity class
  avg.patient.cells <- log1p(AverageExpression(Sobj, verbose = FALSE)$RNA)
  avg.patient.cells.df <- data.frame(avg.patient.cells)
  
  # set the positions for the annotations
  x_axisname_string <- deparse(substitute(x_axis))
  y_axisname_string <- deparse(substitute(y_axis))
  annotation_x <- avg.patient.cells.df[[x_axisname_string]]
  annotation_y <- avg.patient.cells.df[[y_axisname_string]]
  
  
  scatter.plot <- ggplot(avg.patient.cells.df, aes(x = {{ x_axis }}, y = {{ y_axis}} )) +
    geom_point() +
    xlab(paste0(deparse(substitute(x_axis)), " cells, Log2 average gene expression")) + 
    ylab(paste0(deparse(substitute(y_axis)), " cells, Log2 average gene expression")) +
    geom_abline(slope=1, intercept = 0, col="red") +
    geom_abline(slope=1, intercept = 0.6, col="red") +
    geom_abline(slope=1, intercept = -0.6, col="red") +
    annotate("text", x = max(annotation_x) - 1, y = 0.5, label = paste0("Upregulated \n in ", x_axisname_string, " cells"), color = "black") +
    annotate("text", x = 1, y = max(annotation_y) - 0.5, label = paste0("Upregulated \n in ", y_axisname_string,  " cells"), color = "black")
  
  scatter.plot.labels <- LabelPoints(plot = scatter.plot, points = genes_to_label[1:top_genes], repel = TRUE, xnudge=0, ynudge=0)
  
  return(scatter.plot.labels)
}


# function: MA_plot_DEgenes ----------------------------------------------------
# insert a seurat object, the pseudobulk differential expression analysis results table (DE_genes),
# the differentially expressed genes (genes_to_label),
# how many of the differentially expressed genes to be labeled on the plot (top_genes)

MA_plot_DEgenes <- function(Sobj, DE_genes, genes_to_label, top_genes){
  
  geneLogSums <- log2(rowSums(GetAssayData(Sobj, "counts", "RNA")))
  DE_genes$logSum <-geneLogSums[rownames(DE_genes)]
  
  ma.plot <- DE_genes %>%
    ggplot() +
    aes(x = logSum, y = avg_log2FC, col = p_val_adj) +
    geom_point()
  # Add vertical lines for log2FoldChange thresholds
  ma.plot <- ma.plot + geom_hline(yintercept = c(-0.6, 0.6), col = "red") +
    xlab("Log2 Mean Expression") + ylab("Log2 Fold Change") 
  # Add labels to the genes
  ma.plot <- LabelPoints(plot = ma.plot, points = genes_to_label[1:top_genes], repel = TRUE, xnudge=0, ynudge=0)
  
  return(ma.plot)
}


# function: volcano_plot_DEgenes -----------------------------------------------
# insert the pseudobulk differential expression analysis results table (DE_genes),
# the differentially expressed genes (genes_to_label),
# how many of the differentially expressed genes to be labeled on the plot (top_genes),
# and what to be plotted on the y axis: "-log10(p_val_adj)" or "1/p_val_adj"

# WARNING: Seurat rounds very small p values to 0!! 
# Idea on how to approach this:
# Replace all zeros in pval & padjval with 1/2 the minimum value greater than zero from these columns
# DE.genes["p_val"] <- lapply(DE.genes["p_val"], function(x) replace(x, x == 0, val/2 - runif( nrows, 0, val/4)  ))
# DE.genes["p_val_adj"] <- lapply(DE.genes["p_val_adj"], function(x) replace(x, x == 0, val/2 - runif( nrows, 0, val/4)   ))

volcano_plot_DEgenes <- function(DE_genes, genes_to_label, top_genes, y_axis){
  
  if (y_axis == "-log10(p_val_adj)"){
    volcano.plot <- DE_genes %>%
      ggplot() +
      aes(x = avg_log2FC, y = -log10(p_val_adj), col = diffexpressed) +
      scale_y_log10() +
      geom_point()
    
    y_label <- "-log10(p-adjusted-value)"
    yintercept_value <- -log10(0.05)
  }
  
  else if (y_axis == "1/p_val_adj"){
    volcano.plot <- DE_genes %>%
      ggplot() +
      aes(x = avg_log2FC, y= 1 / p_val_adj, col = diffexpressed) +
      scale_y_log10() +
      geom_point()
    
    y_label <- "1/p-adjusted-value"
    yintercept_value <- 1/0.05
  }
  
  # Add vertical lines for log2FoldChange thresholds, and one horizontal line for the p-value threshold
  volcano.plot <- volcano.plot + geom_vline(xintercept = c(-0.6, 0.6), col = "red") +
    geom_hline(yintercept = yintercept_value, col = "red") +
    xlab("Log2 Fold Change") #+ 
    # ylab(paste0("Significance \n ", y_label))  
  
  ## Change point color
  # 1. by default, it is assigned to the categories in an alphabetical order):
  volcano.plot <- volcano.plot + scale_color_manual(values=c("blue", "black", "red"))
  # 2. to automate a bit: create a named vector: the values are the colors to be used, the names are the categories they will be assigned to:
  mycolors <- c("blue", "red", "black")
  names(mycolors) <- c("DOWN", "UP", "NO")
  volcano.plot <- volcano.plot + scale_colour_manual(values = mycolors)
  
  volcano.plot <- LabelPoints(plot = volcano.plot, points = genes_to_label[1:top_genes], repel = TRUE, xnudge=0, ynudge=0)
  
  return(volcano.plot)
}