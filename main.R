# main.R
source("src/data_loader.R")
source("src/qc_normalize.R")
source("src/dim_reduction.R")
source("src/clustering.R")
source("src/enrichment.R")
source("src/plots_visualizer.R")

message("=== Starting Pipeline ===")

# Explicit directory verification and runtime check
results_dir <- "results"
if (!dir.exists(results_dir)) {
  message(paste("Directory '", results_dir, "' not found. Generating path...", sep=""))
  dir.create(results_dir, recursive = TRUE)
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

message("=== Running Heads-Up Graphical Figure Exporter ===")
# Invoke programmatic diagram export
export_pipeline_visualizations(se, results_dir = results_dir)

message("=== Exporting Tabular Data Metrics ===")
markers <- identify_markers(se)
write.csv(markers, file = file.path(results_dir, "cluster_markers.csv"), row.names = FALSE)

# Execute FEA for validation on an arbitrary cluster (e.g., Cluster 0)
fea_results <- run_functional_enrichment(markers, target_cluster = "0")
if (!is.null(fea_results)) {
  write.csv(fea_results, file = file.path(results_dir, "functional_enrichment_cluster0.csv"), row.names = FALSE)
}
message("=== Pipeline Run Completed Successfully ===")