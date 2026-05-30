library(Seurat)

#' Step 1: Programmatic Data Retrieval
download_data <- function(data_dir = "data") {
  if (!dir.exists(data_dir)) {
    dir.create(data_dir, recursive = TRUE)
  }
  
  meta_url <- "https://raw.githubusercontent.com/jdblischak/singleCellSeq/master/data/tung_metadata.tsv"
  counts_url <- "https://raw.githubusercontent.com/jdblischak/singleCellSeq/master/data/molecules.txt"
  
  meta_path <- file.path(data_dir, "tung_metadata.tsv")
  counts_path <- file.path(data_dir, "molecules.txt")
  
  if (!file.exists(meta_path)) {
    message("Downloading Tung et al. metadata...")
    download.file(meta_url, destfile = meta_path)
  }
  if (!file.exists(counts_path)) {
    message("Downloading Tung et al. molecule count matrix...")
    download.file(counts_url, destfile = counts_path)
  }
  return(list(metadata = meta_path, counts = counts_path))
}

#' Steps 2 & 3: Data Ingestion, Formatting, and Identity Alignment
load_and_align_data <- function(file_paths) {
  message("Parsing files into memory...")
  metadata <- read.table(file_paths$metadata, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  counts <- read.table(file_paths$counts, header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
  
  # Align dimensions: ensure intersections match across cell columns and metadata rows
  common_cells <- intersect(rownames(metadata), colnames(counts))
  metadata <- metadata[common_cells, , drop = FALSE]
  counts <- counts[, common_cells]
  
  # Integrity check verification
  if (!all(rownames(metadata) == colnames(counts))) {
    stop("Data alignment failure: Column names do not exactly match metadata row names.")
  }
  
  message("Instantiating Seurat Object...")
  seurat_obj <- CreateSeuratObject(counts = counts, meta.data = metadata, project = "Tung_iPSC")
  return(seurat_obj)
}