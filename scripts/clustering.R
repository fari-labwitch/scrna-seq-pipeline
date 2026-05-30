library(Seurat)

compute_clusters <- function(seurat_obj, dims = 1:30, resolution = 0.5) {
  seurat_obj <- FindNeighbors(seurat_obj, dims = dims, verbose = FALSE)
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution, verbose = FALSE)
  return(seurat_obj)
}

identify_markers <- function(seurat_obj, logfc_threshold = 0.25) {
  return(FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = logfc_threshold, verbose = FALSE))
}