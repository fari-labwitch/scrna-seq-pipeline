library(gprofiler2)
library(dplyr)

#' Run Functional Enrichment Analysis and flatten nested list columns for CSV export
run_functional_enrichment <- function(markers, target_cluster, log2fc_cutoff = 0.5) {
  sig_markers <- markers %>% 
    filter(cluster == target_cluster & avg_log2FC > log2fc_cutoff) %>% 
    pull(gene)
  
  if (length(sig_markers) == 0) {
    warning("No marker features exceeded the selected threshold parameters.")
    return(NULL)
  }
  
  enrichment <- gost(query = sig_markers, organism = "hsapiens", significant = TRUE, sources = "GO:BP")
  
  # Check if any enrichment terms were returned
  if (is.null(enrichment) || is.null(enrichment$result)) {
    warning(paste("No significant GO pathways found for cluster", target_cluster))
    return(NULL)
  }
  
  res <- enrichment$result
  
  # CRITICAL FIX: Convert nested list columns (like 'parents') into flat comma-separated strings
  res <- res %>%
    mutate(across(where(is.list), ~ sapply(., paste, collapse = ", ")))
  
  return(res)
}