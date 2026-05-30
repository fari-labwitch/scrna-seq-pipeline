library(Seurat)

#' Step 1: Programmatic Data Retrieval from the true singleCellSeq repository
download_data <- function(data_dir = "data") {
  if (!dir.exists(data_dir)) {
    dir.create(data_dir, recursive = TRUE)
  }
  
  # Fixed URLs pointing to the exact repository assets
  meta_url <- "https://raw.githubusercontent.com/jdblischak/singleCellSeq/master/data/annotation.txt"
  counts_url <- "https://raw.githubusercontent.com/jdblischak/singleCellSeq/master/data/molecules.txt"
  
  meta_path <- file.path(data_dir, "annotation.txt")
  counts_path <- file.path(data_dir, "molecules.txt")
  
  if (!file.exists(meta_path)) {
    message("Downloading Tung et al. cell annotations...")
    download.file(meta_url, destfile = meta_path, method = "libcurl")
  }
  if (!file.exists(counts_path)) {
    message("Downloading Tung et al. UMI molecule count matrix...")
    download.file(counts_url, destfile = counts_path, method = "libcurl")
  }
  return(list(metadata = meta_path, counts = counts_path))
}

#' Steps 2 & 3: Data Ingestion, Robust Formatting, and Identity Alignment
load_and_align_data <- function(file_paths) {
  message("Parsing files into memory...")
  
  # Read expression count matrix (Setting row.names = 1 captures the Ensembl gene keys)
  counts <- read.table(file_paths$counts, header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
  
  # Read annotations metadata
  metadata <- read.table(file_paths$metadata, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  
  # Robust structural check: If row names don't natively overlap, match using the first tracking column
  if (!any(rownames(metadata) %in% colnames(counts)) && any(metadata[, 1] %in% colnames(counts))) {
    rownames(metadata) <- metadata[, 1]
  }
  
  # Intersect columns and rows to find matching single cells
  common_cells <- intersect(rownames(metadata), colnames(counts))
  
  # Positional index fallback if naming structures vary implicitly
  if (length(common_cells) == 0) {
    if (nrow(metadata) == ncol(counts)) {
      message("Warning: Cell strings did not explicitly cross-match. Syncing indices positionally.")
      rownames(metadata) <- colnames(counts)
      common_cells <- colnames(counts)
    } else {
      stop("Fatal Error: Expression columns cannot be safely aligned with metadata dimensions.")
    }
  }
  
  # Strict alignment slice
  metadata <- metadata[common_cells, , drop = FALSE]
  counts <- counts[, common_cells]
  
  message(paste("Successfully synced alignment for:", length(common_cells), "cells."))
  message("Instantiating Seurat Object...")
  
  seurat_obj <- CreateSeuratObject(counts = as.matrix(counts), meta.data = metadata, project = "Tung_iPSC")
  return(seurat_obj)
}