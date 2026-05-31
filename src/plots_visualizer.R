library(Seurat)
library(ggplot2)

#' Main Orchestration Utility to Save All Analytical Figures
export_pipeline_visualizations <- function(seurat_obj, results_dir = "results") {
  if (!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
  }
  
  # 1. Quality Control Metrics (Violin Distribution)
  # Captures distribution of individual molecules vs raw unique feature genes
  message("-> Exporting QC Metrics Violin Plot...")
  png(filename = file.path(results_dir, "01_qc_metrics_distribution.png"), 
      width = 1000, height = 700, res = 120)
  p1 <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA"), group.by = "replicate") +
    labs(title = "Single-Cell QC Metrics Distribution", subtitle = "Grouped by Technical Replicate Stream")
  print(p1)
  dev.off()
  
  # 2. Gene Counts vs Outlier Filtering Histogram
  message("-> Exporting Unique Feature Ingestion Histogram...")
  png(filename = file.path(results_dir, "02_cell_feature_histogram.png"), 
      width = 800, height = 600, res = 120)
  p2 <- ggplot(seurat_obj[[]], aes(x = nFeature_RNA)) +
    geom_histogram(bins = 50, fill = "#4582ec", color = "black", alpha = 0.7) +
    theme_classic() +
    labs(title = "Distribution of Extracted Genes Per Single Cell",
         x = "Number of Unique Genes (Features)", y = "Frequency (Cell Counts)")
  print(p2)
  dev.off()
  
  # 3. PCA Subspace Dimensional Mapping
  message("-> Exporting PCA Feature Projection...")
  png(filename = file.path(results_dir, "03_pca_dimensionality_reduction.png"), 
      width = 900, height = 700, res = 120)
  p3 <- DimPlot(seurat_obj, reduction = "pca", dims = c(1, 2), group.by = "replicate") +
    labs(title = "Principal Component Analysis (PCA) Coordinates", 
         subtitle = "Evaluating Linear Batch Effects Across Replicates")
  print(p3)
  dev.off()
  
  # 4. Non-Linear UMAP Cluster Map
  message("-> Exporting Graph-Based UMAP Clustering Space...")
  png(filename = file.path(results_dir, "04_umap_cluster_map.png"), 
      width = 900, height = 700, res = 120)
  p4 <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, label.size = 5) +
    labs(title = "Uniform Manifold Approximation & Projection (UMAP)", 
         subtitle = "Identified Unsupervised Cell Populations") +
    theme(plot.title = element_text(face = "bold"))
  print(p4)
  dev.off()
  
  # 5. Top Highly Variable Features Scatter Plot
  if ("vst" %in% names(seurat_obj@assays$RNA@meta.features) || "SCT" %in% names(seurat_obj@assays)) {
    message("-> Exporting Highly Variable Feature Dispersion Map...")
    png(filename = file.path(results_dir, "05_highly_variable_genes.png"), 
        width = 1000, height = 600, res = 120)
    p5 <- VariableFeaturePlot(seurat_obj) + 
      labs(title = "Top Variable Genes Selection Profile", subtitle = "Bioinformatic Variance Stabilization Mode")
    print(p5)
    dev.off()
  }
  
  message("All structural visualization graphs successfully written to destination path.")
}