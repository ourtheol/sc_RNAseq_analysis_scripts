library("enrichR")

library("msigdbr")
library("fgsea")

library("org.Hs.eg.db")
library("AnnotationDbi")
library("clusterProfiler")
library("enrichplot")


# function: gene_set_enrichment ------------------------------------------------
# Gene Set Analysis: perform enrichment using the enrichR package
# input a list of genes, look for their combined function in list of databases
# check package enrichR for available DBs

gene_set_enrichment <- function(genes, databases) {

    enrich_results <- enrichr(genes = genes, databases = databases)
    print(enrich_results)
    
    for (i in length(databases)) {
      
      write.csv(x=enrich_results[[i]], file=paste0(databases[i], "_gene_set_enrichment.csv"))
      
      enrich.plot <- plotEnrich(enrich_results[[i]], showTerms = 25, numChar = 75, y = "Count", orderBy = "P.value")
      
      pdf(paste0(databases[i], "_enrich_plot.pdf"), width = 25, height = 20, bg = "white")
      print(enrich.plot)
      dev.off()
    }
}


# function: Gsea_topPathways ---------------------------------------------------
# check if a particular gene set is more present in the UP-regulated genes, among the DOWN_regulated genes or not differentially regulated
# insert the differential expression analysis results table (DE_resuls) and
# which gene set from msigdbr you want to use (choose among: "CP:KEGG", "CP:WIKIPATHWAYS", "IMMUNESIGDB")

Gsea_topPathways <- function(DE_resuls, gene_set_to_use) { 
 
  # Create a gene rank based on the gene expression fold change
  gene_rank <- setNames(DE_resuls$avg_log2FC, casefold(rownames(DE_resuls),upper = T))

  # Download gene sets
  msigdbgmt <- msigdbr::msigdbr("Homo sapiens")
  msigdbgmt <- as.data.frame(msigdbgmt)
  
  # Subset which gene set you want to use
  msigdbgmt_subset <- msigdbgmt[msigdbgmt$gs_subcat == gene_set_to_use, ]
  gmt <- lapply(unique(msigdbgmt_subset$gs_name), function(x) {
    msigdbgmt_subset[msigdbgmt_subset$gs_name == x, "gene_symbol"]
  })
  names(gmt) <- unique(paste0(msigdbgmt_subset$gs_name, "_", msigdbgmt_subset$gs_exact_source))
  
  # Perform enrichment analysis
  fgseaRes <- fgsea(pathways = gmt, stats = gene_rank, eps = 0.0, minSize = 15, maxSize = 500)
  
  topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=15), pathway]
  topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=15), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

  write.csv(x=topPathways, file="topPathways.csv")
  
  # Nice summary table (shown as a plot)
  pdf("GseaTable.pdf", width = 30, height = 20, bg = "white")
  plotGseaTable(
    pathways = gmt[topPathways],
    stats = gene_rank,
    fgseaRes = fgseaRes,
    gseaParam = 0.5,
    colwidths = c(5, 3, 0.8, 1.2, 1.2),
    render = TRUE)
  dev.off()
 }


# function: enrichGO_plots -----------------------------------------------------
# insert the differential expression analysis results table (DE_resuls) and
# get a pdf with various gene set enrinchment analysis plots 

enrichGO_plots <- function(DE_resuls) { 
  
  # map DE gene names to Entrez IDs --> enrichGO() input
  # drop the non-mapped entries
  DE.genes.plusentrez <- DE_resuls %>%
    mutate(entrez = mapIds(org.Hs.eg.db, keys = DE_resuls$delabel, column = "ENTREZID", keytype = "SYMBOL")) %>%
    filter(!is.na(entrez))
  
  # re-order the dataframe in decreasing order based on the log2FC and then only include columns: ENTREZ and logFC
  DE.entrez.logFC.reordered <- DE.genes.plusentrez[order(DE.genes.plusentrez$avg_log2FC,decreasing = TRUE), c("entrez", "avg_log2FC")]
  
  # create a numeric vector rownames: ENTREZ IDs column entries: log2FC, for cnetplot and heatmap functions in order to visualize the FC
  geneList_PC <- DE.entrez.logFC.reordered$avg_log2FC   
  names(geneList_PC) <- DE.entrez.logFC.reordered$entrez  
  ## omit any NA values
  # geneList_PC <-na.omit(original_geneList_PC)  
  ##  sort again the list in decreasing order (required for clusterProfiler)
  # geneList_PC = sort(geneList_PC, decreasing = TRUE)
  
 
  ego <- enrichGO(DE.entrez.logFC.reordered$entrez, OrgDb = "org.Hs.eg.db", ont = "BP", readable = TRUE)
  
  write.csv(x=ego, file="enrichGO.csv")
  
  barplot <- barplot(ego, showCategory = 25) + theme(text = element_text(size = 8), axis.text.y = element_text(size = 8))&theme(aspect.ratio = 1)

  dotplot <- dotplot(ego, showCategory = 25)+ theme(text = element_text(size = 8),axis.text.y = element_text(size = 8))&theme(aspect.ratio = 1)

  ## remove redundant GO terms
  ego2 <- simplify(ego)
  
  cnet <- cnetplot(ego2, foldChange=geneList_PC, circular=TRUE, colorEdge=TRUE, cex_category=0.7, cex_gene=0.5, cex_label_category=0, cex_label_gene=0.6)+
    theme(legend.title = element_text(size=10), legend.text = element_text(size=10))&theme(aspect.ratio=1)+
    theme(legend.position="bottom", legend.direction="vertical")
  
  heatplot <- heatplot(ego2, foldChange=geneList_PC)+ 
    ggplot2::coord_flip()+
    ggplot2::ylab('Annotations of Biological Processes')+
    ggplot2::xlab('Gene Symbols')+
    ggplot2::theme(axis.text.x=element_text(size=8))+
    ggplot2::theme(axis.text.y=element_text(size=4))+
    ggplot2::ggtitle('Over-representation for Biological Processes')+
    ggplot2::theme(panel.background=element_rect(fill="black",colour="black"))

  
  pdf("enrichGO_plots.pdf")
  print(barplot)
  print(dotplot)
  print(cnet)
  print(heatplot)
  dev.off()
}