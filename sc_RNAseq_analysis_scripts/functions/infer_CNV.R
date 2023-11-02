library("Seurat")
library("infercnv")
# library("dplyr")
# library("tidyr")
library("ggplot2")


# function: inferCNV_input_file_prep -------------------------------------------
# insert the integrated seurat object for the raw_counts_matrix and annotations_file to be created

enrichGO_plots <- function(integrated_samples) {
  # Input 1: Raw counts matrix as tsv file
  write.table(as.matrix(GetAssayData(object = integrated_samples, slot = "counts")),
              'rawcountsmatrix.tsv',
              sep = '\t', row.names = T, col.names = T, quote = F)
  
  # Input 2: Annotation file
  annotation_file <- subset(integrated_samples@meta.data, select = c(orig.ident))
  
  write.table(as.matrix(annotation_file2),
              'annotationfile2.tsv',
              sep = '\t', row.names = T, col.names = F, quote = F)
  
  # Input 3: gene - chromosome positions file
  # wget file: hg38_gencode_v27.txt from https://data.broadinstitute.org/Trinity/CTAT/cnv/
  # hg38.genes <- read.table("hg38_gencode_v27.txt",sep="\t")
}


# function: normalize ----------------------------------------------------------
# Function to normalize column-wise between -1 and 1
# For each cell, gene expression was re-standardized and values were limited as -1 to 1

normalize <- function(x) {
  x <- sweep(x, 2, apply(x, 2, min))
  x <- sweep(x, 2, apply(x, 2, max), "/")
  2*x - 1
}


# function: mutation_burden ----------------------------------------------------
# Calculate the CNV score of each cell, namely mutation burden
# Calculated as the quadratic sum of CNV for each gene.
# input an inferCNV object and the path to the predicted cnv genes (eg: "HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.05.pred_cnv_genes.dat")

mutation_burden <- function(infercnv_obj, path_to_cnv_genes) {
  
  expr.matrix <- infercnv_obj@expr.data
  
  expr.matrix.scaled <- normalize(expr.matrix)
  
  ## We need the quadratic sum of CNVregions, lets collect the genes with CNVs
  CNV.genes.dat<-read.table(path_to_cnv_genes, header=TRUE)
  CNV.genes<-unique(CNV.genes.dat$gene)
  
  ## From the scaled matrix keep only the CNV.genes (rows)
  expr.matrix.scaled.subset <- expr.matrix.scaled[CNV.genes,]
  
  ## The CNV score of each cell, namely mutation burden was calculated as quadratic sum of CNV for each gene.
  # Calculate the quadratic sum of each column
  quadratic_sums <- apply(expr.matrix.scaled.subset, 2, function(col) {
    sum(col^2)
  })
  
  quadratic_sums.df <- as.data.frame(quadratic_sums)
  
  write.table(quadratic_sums.df,
              'quadratic_sums.tsv',
              sep = '\t', row.names = T, col.names = F, quote = F)
  
  return(quadratic_sums.df)
}



# Example of running inferCNV

infercnv_obj = CreateInfercnvObject(raw_counts_matrix="rawcountsmatrix.tsv",
                                    annotations_file="annotationfile.tsv",
                                    delim="\t",
                                    gene_order_file="hg38_gencode_v27.txt",
                                    ref_group_names=c("Healthy_WM0004", "Healthy_WM4072"))

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             ##min number of reads a gene must have on average across all cells not to be filtered out
                             ##could set it to 0 if you have used the raw counts matrix downstream of Seurat
                             out_dir="./",
                             cluster_by_groups=TRUE,
                             denoise=TRUE,
                             HMM=TRUE,
                             # up_to_step=15, ##stop analysis and view intermediate results, 15 -> stop after sub-clustering
                             useRaster=FALSE,
                             # clustering resolution
                             # higher number, higher number of clusters
                             leiden_resolution=0.001,
                             # posterior probability of each CNV to actually not be a CNV
                             # BayesMaxPNormal=0.5(default)
                             # CNV predictions filtered to have more than 0.5 posterior probability of being neutral
                             # stricter filtering: rerun using 0.2 --> fewer CNVs predicted
                             BayesMaxPNormal=0.05,
                             num_threads=16,
                             no_plot=TRUE)

infercnv_obj <- readRDS('run.final.infercnv_obj')

# Plotting
plot_cnv(
  infercnv_obj,
  out_dir = "./",
  title = "inferCNV figure",
  obs_title = "B-cells from patients",
  ref_title = "B-cells from Healthy donors",
  cluster_by_groups = TRUE,
  cluster_references = TRUE,
  plot_chr_scale = FALSE,
  chr_lengths = NULL,
  k_obs_groups = 1,
  contig_cex = 1,
  x.center = mean(infercnv_obj@expr.data),
  x.range = "auto",
  hclust_method = "ward.D",
  custom_color_pal = NULL,
  color_safe_pal = FALSE,
  output_filename = "infercnv",
  output_format = "pdf",
  png_res = 300,
  dynamic_resize = 1,
  ref_contig = NULL,
  write_expr_matrix = FALSE,
  write_phylo = TRUE,
  useRaster = FALSE
)

plot_per_group(
  infercnv_obj,
  on_references = TRUE,
  on_observations = TRUE,
  sample = FALSE,
  n_cells = 1000,
  every_n = NULL,
  above_m = 1000,
  k_obs_groups = 1,
  base_filename = "infercnv_per_group",
  output_format = "png",
  write_expr_matrix = TRUE,
  save_objects = FALSE,
  png_res = 300,
  dynamic_resize = 0,
  out_dir = "./"
)

mutation_burden.df <- mutation_burden(infercnv_obj,
                                      "HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.05.pred_cnv_genes.dat")
# continue by groupping the cells and plot violin plots of the infer CNV scores, determine statistically significant differences