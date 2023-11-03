library("Seurat")
# Bone Marrow Dataset is available through SeuratData in library: SeuratData
library("SeuratData")
# # Install the reference data
# InstallData("bmcite")


# Mapping to (healthy) human bone marrow cells reference data


# function: bm_reference_prep --------------------------------------------------
# prepares the bone marrow reference dataset

bm_reference_prep <- function() {

  # Load reference data
  bm <- LoadData(ds = "bmcite")
  
  # The reference dataset contains a WNN graph, reflecting a weighted combination of the RNA and protein data in this CITE-seq experiment.
  # We can compute a UMAP visualization based on this graph. We set return.model = TRUE, which will enable us to project query datasets onto this visualization.
  bm <- RunUMAP(bm, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", return.model = TRUE)
  
  # pdf("mapping_reference_bm_dataset_umap.pdf", width = 10, height = 10)
  # DimPlot(bm, group.by = "celltype.l2", reduction = "wnn.umap")&theme(aspect.ratio=1)
  # dev.off()
  
  # The sPCA calculation is performed once, and then can be rapidly projected onto each query dataset.
  bm <- ScaleData(bm, assay = 'RNA')
  
  bm <- RunSPCA(bm, assay = 'RNA', graph = 'wsnn')
  
  bm <- FindNeighbors(
    object = bm,
    reduction = "spca",
    dims = 1:50,
    graph.name = "spca.annoy.neighbors",
    k.param = 50,
    cache.index = TRUE,
    return.neighbor = TRUE,
    l2.norm = TRUE
  )
  
  # bm[["spca.annoy.neighbors"]]
  SaveAnnoyIndex(object = bm[["spca.annoy.neighbors"]], file = "./reftmp.idx")
  # bm[["spca.annoy.neighbors"]] <- LoadAnnoyIndex(object = bm[["spca.annoy.neighbors"]], file = "./reftmp.idx")
  
  return(bm)
}


# function: bm_mapping ---------------------------------------------------------
# input the reference (bone marrow) and a query sample (seurat object) to be mapped on the reference.

# WARNING: Normalize the query in the same manner as the reference. 
# Above, the reference was normalized using log-normalization via NormalizeData(). 
# If the reference had been normalized using SCTransform(), the query must be normalized with SCTransform() as well.

# returns the sample (seurat object) with added meta.data column: predicted.celltype 

bm_mapping <- function(bm, sample) {
  
  # Find anchors between the donor query dataset and the multimodal reference
  anchors <- FindTransferAnchors(
    reference = bm,
    query = sample,
    k.filter = NA,
    reference.reduction = "spca",
    reference.neighbors = "spca.annoy.neighbors",
    dims = 1:50
  )
  
  # Map the dataset
  sample <- MapQuery(
    anchorset = anchors,
    query = sample,
    reference = bm,
    refdata = list(
      celltype = "celltype.l2",
      predicted_ADT = "ADT"),
    reference.reduction = "spca",
    reduction.model = "wnn.umap"
  )
  
  return(sample)
}
