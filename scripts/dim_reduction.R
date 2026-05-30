library(Seurat)

run_dimensionality_reduction <- function(seurat_obj, n_features = 2000, npcs = 30) {
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = n_features, verbose = FALSE)
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
  seurat_obj <- RunPCA(seurat_obj, npcs = npcs, verbose = FALSE)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:npcs, verbose = FALSE)
  return(seurat_obj)
}