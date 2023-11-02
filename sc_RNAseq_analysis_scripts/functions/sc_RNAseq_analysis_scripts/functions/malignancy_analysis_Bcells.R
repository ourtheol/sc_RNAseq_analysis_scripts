library("dplyr")
library("tidyr")
library("Seurat")
library("ggplot2")
library("cowplot")
library("purrr")
library("ggrepel")
library("tidyverse")
library("RColorBrewer")
library("Seurat")
library("grid")


# function: get.data -----------------------------------------------------------
# insert a seurat object, and
# return a data frame with clustering results, dimensionality reduction values, selected genes expression

get.data <- function(Sobj, genes) {
  
  # Extract Meta data
  # insert column "Cluster" with the current cell identities (usually the column: seurat_clusters)
  df.meta <- data.frame(Cluster = Sobj@active.ident, Sobj@meta.data)
  
  # Extract tSNE coordinates if available
  # insert columns "tSNE_1" and "tSNE_2"
  if(!is.null(Sobj@reductions$tsne))
  {
    df.meta <- cbind(df.meta, Sobj@reductions$tsne@cell.embeddings)
  }
  
  # Extract UMAP coordinates if available
  # insert columns "UMAP_1" and "UMAP_2"
  if(!is.null(Sobj@reductions$umap))
  {
    df.meta <- cbind(df.meta, Sobj@reductions$umap@cell.embeddings)
  }
  
  # Extract gene expression data for expressed genes
  # insert columns of expression of selected genes (IGLC2,IGKC)
  genes.red <- rownames(Sobj@assays$RNA@data)[rownames(Sobj@assays$RNA@data) %in% genes]
  
  dftotal <-  cbind(df.meta, FetchData(Sobj, genes.red))

  return(dftotal)
}


# function: KL_ratio calculation -----------------------------------------------
# import the previously created df with the columns extracted from a seurat object, and
# calculate the IGKC-fraction of each cell

KL_ratio <- function(x) {
  
  df <- get.data(x, genes=c("IGLC2", "IGKC"))
  
  if(!"IGKC" %in% colnames(df)) { df <- mutate(df, IGKC=0)}
  if(!"IGLC2" %in% colnames(df)) { df <- mutate(df, IGLC2=0)}
  
  # To avoid overplotting, ratios of 1 are jittered with a negative value,
  # whereas values of 0 are jittered with a positive value.
  
  # insert Ratio.jittered column, the IGKC-fraction of each cell
  df <- df %>% mutate(Ratio.jittered = ifelse(IGKC/(IGLC2+IGKC) == 1,
                                              IGKC/(IGLC2+IGKC) - runif( n(), 0, 0.2 ),
                                              ifelse(IGKC/(IGLC2+IGKC) == 0,
                                                     IGKC/(IGLC2+IGKC) + runif( n(), 0, 0.2 ), IGKC/(IGLC2+IGKC))))
  
  df <- df %>% mutate(kappa.lambda.ratio = IGKC/(IGLC2+IGKC))
  
  return(df)
}


# function: Plot.KL, plots from Roider et al., 2020 ----------------------------
# https://github.com/DietrichLab/scLymphomaExplorer/blob/master/Rmd/scLN.Analysis.Rmd
# Identification of light chain restricted malignant cell populations:
# Light chain expression in non-malignant cells results in a plot which is blue/red scattered. 
# Single-colored clusters indicate light chain restricted, therefore malignant cells. 

# insert a data frame with clustering results, dimensionality reduction values, IGKC and IGLC2 expression and the IGKC-fraction of each cell
# and select the clustering results column to be plotted

plot.KL <- function(df, column){
  
  # UMAP dimentionality reduction plot, colors according to each cell's cluster 
  p1 <- ggplot()+
    geom_point(data = df,
               aes(x = UMAP_1, y = UMAP_2, color = column), 
               size=0.5, alpha=0.75)+
    theme_bw()+
    theme(aspect.ratio = 1)+
    guides(color = guide_legend(override.aes = list(size = 2))) +
    xlab("UMAP1")+ylab("UMAP2")
  
  # Scater plot: place cells according to their expression values along the IGKC and IGLC2 axis
  p2 <- ggplot(df)+
    geom_rect(aes(xmin = -Inf, xmax = +Inf, ymin = -0.30, max = 0.30), fill = "grey95", alpha = 0.05)+
    geom_rect(aes(xmin = -0.30, xmax = 0.30, ymin = -Inf, max = Inf), fill = "grey95", alpha = 0.05)+
    geom_jitter(size=0.75, aes(x = IGKC, y = IGLC2, color = column), height = 0.25, width = 0.25)+   
  ylim(c(-0.5, 7.5))+
    xlim(c(-0.5, 7.5))+
    theme_bw()+
    guides(color = guide_legend(override.aes = list(size = 2))) +
    theme(axis.title = element_text(face = "italic"),
          aspect.ratio = 1)
  
  # UMAP dimentionality reduction plot, colors according to each cell's IGKC-fraction value
  p3 <- ggplot()+
    geom_point(data=df,
               aes(x=UMAP_1, y=UMAP_2, color=(IGKC+0.001)/(IGKC+IGLC2+0.001)),
               size=0.75, alpha=0.75)+
    scale_color_gradient2(low = "#fc8d59", mid = "#ffffbf", limits=c(0, 1),
                          midpoint = 0.5, high = "#91bfdb", name = "K/L ratio", guide = "colourbar")+
    theme_bw()+
    theme(aspect.ratio = 1)+
    xlab("UMAP1")+ylab("UMAP2")
  
  # IGKC-fraction for each cell along y axis, cells grouped in clusters 
  p4 <- ggplot()+
    geom_jitter(data=df, size=0.75, aes(x=column, y=Ratio.jittered,    
                                        color=Ratio.jittered), width = 0.25)+
    scale_color_gradient2(low = "#f46d43", mid = "#ffffbf", high = "steelblue2", midpoint = 0.5,
                          na.value = "#fc8d59", guide = "none")+
    ylim(c(0, 1))+
    ylab(expression(italic(IGKC/(IGLC2+IGKC))))+
    geom_hline(yintercept = 0.5, linetype="dashed")+
    theme_bw()
  theme(aspect.ratio = 1,
        axis.text.y = element_text(size=11),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=12, angle=45, hjust = 1))
  
  
  plot_grid(p1, p3, p2, p4, align = "v", axis = "left", rel_widths = c(1, 1))
}


# function: determine_healthy_neoplastic ---------------------------------------
# insert a data frame with clustering results, dimensionality reduction values, IGKC and IGLC2 expression and the IGKC-fraction of each cell
# and select the clustering results column name (seurat_clusters, RNA_snn_res.0.4, etc) to be plotted 

determine_healthy_neoplastic <- function(ratio, clustering_colname, healthy_lower_level, healthy_upper_level, file_prefix){
  
  # calculate the IGKC-fraction for each cell, then
  # group by selected clustering (seurat_clusters, RNA_snn_res.0.4, etc) and calculate the kappa/lambda ratio for each cluster
  message("Calculating the kappa / lambda ratio for each cluster ...")
  ratio.percluster <- ratio %>%
    # assign kappa or lambda (or healthy) to each cell
    dplyr::mutate(Bcell.type = case_when(kappa.lambda.ratio > 0.5 ~ "kappa",
                                         kappa.lambda.ratio < 0.5 ~ "lambda",
                                         kappa.lambda.ratio == 0.5 ~ "healthy")) %>%
    # count the number of kappa, lambda, healthy cells in each cluster
    dplyr::select({{ clustering_colname }}, Bcell.type) %>% 
    dplyr::group_by({{ clustering_colname }}, Bcell.type) %>%  
    dplyr::count() %>%
    # calculate the percentages of kappa, lambda, healthy in each cluster
    dplyr::group_by({{ clustering_colname }}) %>% 
    dplyr::mutate(per = 100 * n / sum(n)) %>%
    # rotate the matrix to have separate columns of kappa, lambda, healthy
    dplyr::select({{ clustering_colname }},	Bcell.type,	per) %>%  
    pivot_wider(names_from = Bcell.type, values_from = per) %>%
    # replace NA with 1 to facilitate the divisions
    replace(is.na(.), 1) %>%
    # calculate the kappa / lambda ratio for each cluster
    dplyr::mutate(cluster.kappa.lambda.ratio = kappa / lambda) %>%
    dplyr::mutate(cluster.char = if_else(cluster.kappa.lambda.ratio < healthy_lower_level | cluster.kappa.lambda.ratio > healthy_upper_level, "neoplastic", "healthy"))
  write.csv(x=ratio.percluster, file=paste0(file_prefix,  "_kappa_lambda_ratio_percluster.csv")) 
  
  ratio.percell <- left_join(ratio, ratio.percluster, by = deparse(substitute(clustering_colname)))  %>%  #clustering_colname as string
    dplyr::select(orig.ident, {{ clustering_colname }}, kappa.lambda.ratio, cluster.char) 
  write.csv(x=ratio.percell, file=paste0( file_prefix,  "_kappa_lambda_ratio_percell.csv")) 
  
  return(ratio.percell)
}


# function: plot_sample_purity -------------------------------------------------
# insert a data frame with inf per cell: sampleID, clustering, IGKC-fraction, neoplastic/healthy
# plots of the estimated purity of each sample with 95% confidence intervals

plot_sample_purity <- function(meta.data.df){
  
  # Bar plot of purity: count of neoplastic and healthy cells in each sample
  barplot_malignancy <- ggplot(meta.data.df,aes(x=orig.ident, fill=cluster.char))+
    geom_bar()+
    theme(axis.text.x=element_text(angle=45, hjust=1))
  
  # Tumor purity dot plot with confidence intervals
  purity.df <- meta.data.df  %>%
    #count the number of neoplastic cells in each sample
    dplyr::group_by(orig.ident, cluster.char) %>%
    dplyr::count() %>%
    #transpose to have columns: orig.ident, healthy, neoplastic
    pivot_wider(names_from = cluster.char, values_from = n) %>%
    #replace NA with 0
    replace(is.na(.), 0) %>%
    #count the number of all cells in each sample
    dplyr::mutate(total.cells = healthy + neoplastic) %>%
    dplyr::mutate(purity = neoplastic / total.cells) %>%
    dplyr::mutate(margin = qnorm(0.975) * sqrt(purity * (1 - purity) / total.cells)) %>%
    dplyr::mutate(lower = purity - margin) %>%
    dplyr::mutate(upper = purity + margin)
  
  purity_dotplot <- ggplot(data = purity.df, aes(x = orig.ident, y = purity, ymin = lower, ymax = upper)) +
    geom_point(position = position_dodge(width = 0.2)) +
    geom_errorbar(position = position_dodge(width = 0.2), width = 0.1) +
    theme(axis.text.x=element_text(angle = 45, hjust = 1))
  
  # arrange the two plots into one column
  plot_grid(
    barplot_malignancy, purity_dotplot,
    labels = "AUTO", ncol = 1
  )
}











