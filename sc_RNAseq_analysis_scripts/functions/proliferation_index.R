library("dplyr")




# function: proliferation_index ------------------------------------------------
# Calculation of Proliferation Index per cluster
# input a seurat object and the name of the clustering column the indices will be calculated on

proliferation_index <- function(Sobj, clusters) {

  metadata.subset <- Sobj@meta.data %>% 
    dplyr::select({{ clusters }}, Phase) 
  
  # group per cluster , calculate the sum of cells in S, G2M, G1, calculate the indeces:
  
  # Group by count using dplyr
  phase.count.percluster <- metadata.subset %>%
    group_by( {{ clusters }}, Phase) %>%
    summarise(total_count=n(), .groups = 'drop') %>%
    as.data.frame()
  
  ## From long data to wide data
  data_wide <- spread(phase.count.percluster, Phase, total_count) %>%
    # Replace NAs on multiple columns with 0
    mutate(G1 = coalesce(G1, 0),
           G2M = coalesce(G2M, 0),
           S = coalesce(S, 0)) %>%
    mutate( Mitotic_index = G2M  / ( G1 + S + G2M) ) %>%
    mutate( PI = (G2M + S)  / ( G1 + S + G2M) )
  
  write.csv(x=data_wide, file = paste0(Sobj$orig.ident[[1]], "_prolifIndices.csv"))
  
  return(data_wide)
}


# function: plot_proliferation_index ------------------------------------------------
# input the output table of the proliferation_index function and the index to plot

plot_proliferation_index <- function(proliferation_index_out, clusters, index) {
  
  p <- proliferation_index_out %>%
    ggplot( aes(x = {{ clusters }}, y = {{ index }}, group = 1)) +
    geom_line() +
    ## Annotate each cluster by inserting a column: cluster_char, and color differently
    # geom_point(aes(colour = factor(cluster_char)), size = 4) +
    # labs(color='cluster \n characterization') +
    geom_path() +
    ylab("Proliferation Index") +
    xlab("patient clusters") +
    theme(legend.position="none")
  
  return(p)
}