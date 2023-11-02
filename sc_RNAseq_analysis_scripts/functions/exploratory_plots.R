library(Seurat)
library(ggplot2)
library(cowplot)


# function: gene_counts_per_cell_plots -----------------------------------------
# input a seurat object and visualize the number of genes per cell per sample in histogram and boxplot

gene_counts_per_cell_plots <- function(integrated_cells) { 
  
  ## Gene counts per cell
  # write.table(integrated_cells@assays[["RNA"]]@counts, file='Gene_Count_per_Cell.tsv', quote=FALSE, sep='\t', col.names = TRUE)
  ## number of genes per cell: integrated_cells@meta.data$nFeature_RNA
  
  # Visualize the distribution of genes detected per cell via histogram
  gene_distribution_histogram <- integrated_cells@meta.data %>%
    ggplot(aes(color = orig.ident, x = nFeature_RNA, fill = orig.ident)) + 
    geom_density(alpha = 0.2) + 
    theme_classic() +
    scale_x_log10() + 
    geom_vline(xintercept = 300) +
    xlab("number of genes detected per cell") +
    theme(legend.title = element_blank())&theme(aspect.ratio = 1)
  
  # Visualize the distribution of genes detected per cell via boxplot
  # y axis: number of genes
  gene_distribution_boxplot <- integrated_cells@meta.data %>% 
    ggplot(aes(x = orig.ident, y = nFeature_RNA, fill = orig.ident)) + 
    geom_boxplot() + 
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    theme(legend.title=element_blank()) +
    ylab("nGenes")&theme(aspect.ratio = 1)
  
  # Visualize the distribution of genes detected per cell via boxplot
  # y axis: log10(number of genes)
  gene_distribution_boxplot2 <- integrated_cells@meta.data %>% 
    ggplot(aes(x = orig.ident, y = log10(nFeature_RNA), fill = orig.ident)) + 
    geom_boxplot() + 
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    theme(legend.title=element_blank()) +
    ylab("log10(nGenes)")&theme(aspect.ratio = 1)
  
  
  plot_grid(gene_distribution_histogram, gene_distribution_boxplot, gene_distribution_boxplot2, align = "h", nrow = 1)
  }


# function: cell_counts_per_sample_plots -----------------------------------------
# input a seurat object and visualize the number of cells per sample
# option to group and color cells according to a meta.data column

gene_counts_per_cell_plots <- function(integrated_cells, condition) { 
  
  # Bar plot of the number of cells in each sample
  barplot.cells <- ggplot(integrated_cells@meta.data,
                          aes(x = integrated_cells@meta.data$orig.ident, fill = orig.ident))+
    geom_bar()+
    xlab("Sample IDs")+
    ylab("Number of cells")+
    theme(legend.title=element_blank()) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    theme(legend.position = "none")&theme(aspect.ratio = 1)
  
  # Bar plot of the number of cells in each sample colored by a condition
  barplot.condition <- ggplot(integrated_cells@meta.data,
                              aes(x = integrated_cells@meta.data$orig.ident, fill = {{ condition }}))+  
  geom_bar()+
    xlab("Sample IDs")+
    ylab("Number of cells")+
    theme(legend.title = element_blank()) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))&theme(aspect.ratio=1)
  
  # Bar plot of the number of cells in each sample grouped by a condition
  barplot.condition2 <-integrated_cells@meta.data %>%
    group_by({{ condition }}) %>%
    summarise(n=n()) %>%
    ggplot(aes(x = {{ condition }}, y = n, fill = {{ condition }})) + 
    geom_bar(stat = "identity") +
    theme(legend.title = element_blank()) +
    xlab("")+
    ylab("Number of cells")&theme(aspect.ratio = 1)
  
  # Bar plot of the number of cells in each sample grouped by a condition and colored according to each sample
  barplot.condition3 <-integrated_cells@meta.data %>%
    ggplot(aes(x = {{ condition }}, fill = orig.ident)) +
    geom_bar()+
    theme(legend.title = element_blank()) +
    xlab("")+
    ylab("Number of cells")&theme(aspect.ratio = 1)
  
  plot_grid(barplot.cells,
            barplot.condition,
            barplot.condition2,
            barplot.condition3,
            align = "v",
            axis = "left",
            rel_widths = c(1, 1))
  }



## Example of renaming clusters according to cell types
## Set to the clusters from Louvain res 0.4
# Idents(object = integrated_b_cells) <- "RNA_snn_res.0.4"
# levels(integrated_b_cells)
# integrated_b_cells <- RenameIdents(integrated_b_cells,
#                                    "0" = "Non-switched memory B cells",
#                                    "1" = "Non-switched memory B cells",
#                                    "7" = "Non-switched memory B cells",
#                                    "3" = "Switched memory B cells",
#                                    "4" = "Switched memory B cells",
#                                    "9" = "Switched memory B cells",
#                                    "2" = "Naive B cells",
#                                    "8" = "pro-B",
#                                    "6" = "pre-B",
#                                    "5" = "Plasmablasts")
# levels(integrated_b_cells)