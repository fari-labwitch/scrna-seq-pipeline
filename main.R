# main.R
source("scripts/data_loader.R")
source("scripts/qc_normalize.R")
source("scripts/dim_reduction.R")
source("scripts/clustering.R")
source("scripts/enrichment.R")

message("=== Starting Pipeline ===")

# Explicit directory verification and runtime check
results_dir <- "results"
if (!dir.exists(results_dir)) {
  message(paste("Directory '", results_dir, "' not found. Generating path...", sep=""))
  dir.create(results_dir, recursive = TRUE)
} else {
  message(paste("Verified destination directory '", results_dir, "' exists.", sep=""))
}

# Run data ingestion utilities
data_paths <- download_data()
se <- load_and_align_data(data_paths)

message("=== Running Processing and Normalization ===")
se <- apply_qc_filter(se)
se <- normalize_expression(se)

message("=== Computing Dim-Reduction and Clusters ===")
se <- run_dimensionality_reduction(se)
se <- compute_clusters(se, resolution = 0.5)

message("=== Exporting Results ===")
markers <- identify_markers(se)
write.csv(markers, "data/cluster_markers.csv", row.names = FALSE)

# Execute FEA for validation on an arbitrary cluster (e.g., Cluster 0)
fea_results <- run_functional_enrichment(markers, target_cluster = "0")
if (!is.null(fea_results)) {
  write.csv(fea_results, file = file.path(results_dir, "functional_enrichment_cluster0.csv"), row.names = FALSE)
}
message("=== Pipeline Run Completed Successfully ===")