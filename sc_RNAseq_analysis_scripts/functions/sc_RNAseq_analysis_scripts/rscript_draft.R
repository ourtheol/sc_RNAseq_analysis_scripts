setwd("/home/rania/Documents/sc_RNAanalysis_scripts/")
data_path <- ""  ##sometimes it is used as the output path too, correct this
out_path <-  ""


source("functions/scrnaseq_initialanalysis.R")

# Merge samples -----
## merge multiple lists of Seurat objects
samples_merged <- MergeSeurat(samples_list)

saveRDS(samples_merged, "./samples_merged.rds")
samples_merged <- readRDS("samples_merged.rds")


# Filter samples ----

print("Filtering is starting!")

samples_filtered <- FilterSeurat(samples_merged, min_genes =  200, max_genes =  2500, max_mito_genes =  15, min_non_zero_cells_per_gene = 10)

saveRDS(samples_filtered, "./samples_merged_filtered.rds")

print(table(samples_filtered$orig.ident))
samples_filtered <- readRDS("samples_merged_filtered.rds")

# Plot quality control before and after filtering

PlotQC(seurat_merged = samples_merged, seurat_filtered = samples_filtered, "./persample")



# Separate sample analysis, differential expression and plotting ----

samples_analyzed <- SeparateSeuratAnalysis(samples_filtered, 2000, out_path)

saveRDS(samples_analyzed, "./persample_analyzed.rds")
samples_analyzed <- readRDS("persample_analyzed.rds")

samples_degenes_separate <- DifferentiallyExpressedGenes(seurat_analyzed = samples_analyzed, clustering = "seurat_clusters")

saveRDS(samples_degenes_separate, "./persample_degenes_separate.rds")
samples_degenes_separate <- readRDS("persample_degenes_separate.rds")

PrintPlots(samples_degenes_separate, out_path)





source("functions/malignancy_analysis_Bcells.R")

# List of seurat objects that have been processed
samples_analyzed

# Keep data from each sample for each cell (sampleID, clustering, IGKC-fraction, neoplastic/healthy)
meta.data.df <- data.frame()   

for (seurat in 1:length(samples_analyzed)) {
  
  tryCatch({
    
    sampleID <- samples_analyzed[[seurat]]$orig.ident[[1]]
    message("CLONALITY analysis for sample: ")
    message(sampleID)
    
    ratio <- as.data.frame(KL_ratio(samples_analyzed[[seurat]]))
    
    clonality_analysis.leiden.clustering <- plot.KL(ratio, ratio$seurat_clusters)   # or RNA_snn_res.0.4, etc
    
    ratio.percell <- determine_healthy_neoplastic(ratio=ratio, 
                                 clustering_colname = seurat_clusters,   # or RNA_snn_res.0.4, etc 
                                 healthy_lower_level=0, healthy_upper_level=0.01, 
                                 file_prefix = sampleID)
    
    
    meta.data.df <- rbind(meta.data.df, ratio.percell)
    
    # Add the characterization neoplastic or healthy to the Seurat object
    samples_analyzed[[seurat]] <- AddMetaData(object=samples_analyzed[[seurat]], col.name = "neoplastic.or.healthy", metadata = ratio.percell$cluster.char)
    
    message("Printing plots... ")
    message(sampleID)
    pdf(paste0(out_path, sampleID, "_clonality_analysis.pdf"), width = 30, height = 15)
    print(clonality_analysis.leiden.clustering)
    dev.off()
    
    remove(ratio)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

}

# Plot the estimated purity of each sample with 95% confidence intervals
plot_sample_purity(meta.data.df)

saveRDS(data, file = "data.rds")




# make this into a function!! 
# ################################################################################
# # MULTIMODAL REFERENCE MAPPING
# # Mapping human bone marrow cells
# # Bone Marrow Dataset is available through SeuratData in library: SeuratData
# # Load reference data
# # InstallData("bmcite")
# bm <- LoadData(ds = "bmcite")
# # The reference dataset contains a WNN graph, reflecting a weighted combination of the RNA and protein data in this CITE-seq experiment.
# # We can compute a UMAP visualization based on this graph. We set return.model = TRUE, which will enable us to project query datasets onto this visualization.
# bm <- RunUMAP(bm, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", return.model = TRUE)
# # pdf(paste0(out_path, "mapping_reference_bm_dataset_umap.pdf"), width = 10, height = 10)
# # DimPlot(bm, group.by = "celltype.l2", reduction = "wnn.umap")&theme(aspect.ratio=1)
# # dev.off()
# # The sPCA calculation is performed once, and then can be rapidly projected onto each query dataset.
# bm <- ScaleData(bm, assay = 'RNA')
# bm <- RunSPCA(bm, assay = 'RNA', graph = 'wsnn')
# 
# bm <- FindNeighbors(
#   object = bm,
#   reduction = "spca",
#   dims = 1:50,
#   graph.name = "spca.annoy.neighbors",
#   k.param = 50,
#   cache.index = TRUE,
#   return.neighbor = TRUE,
#   l2.norm = TRUE
# )
# 
# # bm[["spca.annoy.neighbors"]]
# SaveAnnoyIndex(object = bm[["spca.annoy.neighbors"]], file = "./reftmp.idx")
# # bm[["spca.annoy.neighbors"]] <- LoadAnnoyIndex(object = bm[["spca.annoy.neighbors"]], file = "./reftmp.idx")
# 
# # Find anchors between each donor query dataset and the multimodal reference
# anchors <- list()
# for (i in 1:length(samples_analyzed)) {
#   anchors[[i]] <- FindTransferAnchors(
#     reference = bm,
#     query = samples_analyzed[[i]],
#     k.filter = NA,
#     reference.reduction = "spca",
#     reference.neighbors = "spca.annoy.neighbors",
#     dims = 1:50
#   )
# }
# # Individually map each of the datasets
# for (i in 1:length(samples_analyzed)) {
#   samples_analyzed[[i]] <- MapQuery(
#     anchorset = anchors[[i]],
#     query = samples_analyzed[[i]],
#     reference = bm,
#     refdata = list(
#       celltype = "celltype.l2",
#       predicted_ADT = "ADT"),
#     reference.reduction = "spca",
#     reduction.model = "wnn.umap"
#   )
# }
# 
# saveRDS(samples_analyzed, file = "datapersample_bm_mapped.rds")




###############################

#  integration part is MISSING!!!!

###############################



source("functions/pseudobulk_Differential_Expression.R")

# data is a seurat object after integration

# Remove from the analysis the imunoglobulin genes
immunoglobulin.LC.genes <- grep("^IGL",rownames(data@assays$RNA@counts),value = TRUE)
immunoglobulin.HC.genes <- grep("^IGH",rownames(data@assays$RNA@counts),value = TRUE)
immunoglobulin.K.genes <- grep("^IGK",rownames(data@assays$RNA@counts),value = TRUE)
immunoglobulin.genes <- c(immunoglobulin.LC.genes,immunoglobulin.HC.genes,immunoglobulin.K.genes)


patients.DE.genes.to.print <- pseudobulk_differential_expession(Sobj = data, 
                                                                condition = "neoplastic.or.healthy", 
                                                                remove_genes = c(), 
                                                                ident_1 = "neoplastic", ident_2 = "healthy")  ## or immunoglobulin.genes

## Get a vector of the differentially expressed genes:
patients.DE.genes.list <- patients.DE.genes.to.print$delabel
print(patients.DE.genes.list)

sc <- scatter_plot_DEgenes(Sobj = data, 
                           condition = "neoplastic.or.healthy", 
                           x_axis = healthy, 
                           y_axis = neoplastic, 
                           genes_to_label = patients.DE.genes.list, 
                           top_genes = 20)


ma <- MA_plot_DEgenes(Sobj = data, patients.DE.genes.to.print,
                      genes_to_label = patients.DE.genes.list, 
                      top_genes = 20)

g <- volcano_plot_DEgenes(patients.DE.genes.to.print,
                          genes_to_label = patients.DE.genes.list, 
                          top_genes = 20,
                          y_axis = "1/p_val_adj")



source("functions/enrichment_analysis.R")
patients.DE.genes.list
head(patients.DE.genes.to.print)
                                                                           
gene_set_enrichment(patients.DE.genes.list, c("GO_Biological_Process_2021","WikiPathway_2021_Human"))
Gsea_topPathways(patients.DE.genes.to.print, "IMMUNESIGDB")

enrichGO_plots(patients.DE.genes.to.print)



source("functions/exploratory_plots.R")
p <- gene_counts_per_cell_plots(MYD88wt_WM0002)
h <- gene_counts_per_cell_plots(MYD88wt_WM0002, condition)


# calculate the proliferation indices on a list of seurat objects
for (seurat in 1:length(data)) {
  
  message(data[[seurat]]@meta.data$orig.ident[[1]])
  
  pro <- proliferation_index(MYD88wt_WM0002, seurat_clusters)
  p <- plot_proliferation_index(pro, seurat_clusters, Mitotic_index)
}

