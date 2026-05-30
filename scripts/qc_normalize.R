library(Seurat)

apply_qc_filter <- function(seurat_obj, min_features = 200, max_features = 10000, max_mito_pct = 5) {
  # Compute mitochondrial gene expression percentages (Human HGNC symbols pattern)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  # Filter outliers representing low-quality or dying single cells
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > min_features & 
                         nFeature_RNA < max_features & 
                         percent.mt < max_mito_pct)
  return(seurat_obj)
}

normalize_expression <- function(seurat_obj, method = "LogNormalize", scale_factor = 10000) {
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = method, scale.factor = scale_factor, verbose = FALSE)
  return(seurat_obj)
}