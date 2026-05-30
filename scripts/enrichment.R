library(gprofiler2)
library(dplyr)

run_functional_enrichment <- function(markers, target_cluster, log2fc_cutoff = 0.5) {
  sig_markers <- markers %>% 
    filter(cluster == target_cluster & avg_log2FC > log2fc_cutoff) %>% 
    pull(gene)
  
  if (length(sig_markers) == 0) {
    warning("No marker features exceeded the selected threshold parameters.")
    return(NULL)
  }
  
  enrichment <- gost(query = sig_markers, organism = "hsapiens", significant = TRUE, sources = "GO:BP")
  return(enrichment$result)
}