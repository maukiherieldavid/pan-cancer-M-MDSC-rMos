# Shared setup for the pan-cancer M-MDSC / rMos analysis.

options(stringsAsFactors = FALSE)
set.seed(1234)

required_packages <- c(
  "Seurat", "dplyr", "tidyr", "tibble", "stringr", "readr",
  "ggplot2", "Matrix", "patchwork"
)

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Install required R packages before running the analysis: ",
       paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(readr)
  library(ggplot2)
  library(Matrix)
  library(patchwork)
})

if (!file.exists(file.path("config", "paths.R"))) {
  stop("Copy config/paths.example.R to config/paths.R and set the local input paths.")
}
source(file.path("config", "paths.R"))

OBJECT_DIR <- file.path(RESULTS_DIR, "objects")
TABLE_DIR <- file.path(RESULTS_DIR, "tables")
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OBJECT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

signature_genes <- c("VCAN", "APOBEC3A", "CD300E", "EREG", "FCN1")
monocytic_context_genes <- c("CD14", "FCN1")
mdsc_context_genes <- c("HLA-DRA", "HLA-DRB1", "TGFB1")
rmo_genes <- c("REL", "CCR2", "NFKB1", "IL1B")

rename_idents <- function(object, mapping) {
  mapping_names <- names(mapping)
  mapping <- as.character(mapping)
  names(mapping) <- mapping_names
  missing_levels <- setdiff(names(mapping), levels(Idents(object)))
  if (length(missing_levels) > 0) {
    stop("Cluster IDs missing from object: ", paste(missing_levels, collapse = ", "))
  }
  do.call(RenameIdents, c(list(object = object), as.list(mapping)))
}

add_qc_metrics <- function(object) {
  object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")
  object <- PercentageFeatureSet(object, pattern = "^HBA|^HBB", col.name = "pHB")
  object <- PercentageFeatureSet(object, pattern = "^RPS|^RPL", col.name = "pRP")
  object
}

normalize_pca <- function(object, nfeatures = 2000, npcs = 50) {
  object <- NormalizeData(object)
  object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = nfeatures)
  object <- ScaleData(object)
  object <- RunPCA(object, features = VariableFeatures(object), npcs = npcs)
  object
}

rank_positive_markers <- function(marker_table, n = 25) {
  marker_table %>%
    rownames_to_column("gene") %>%
    filter(!grepl("^RPS|^RPL", gene)) %>%
    arrange(desc(avg_log2FC)) %>%
    slice_head(n = n)
}

save_table <- function(x, filename) {
  write.csv(x, file.path(TABLE_DIR, filename), row.names = FALSE)
}

save_plot <- function(plot, filename, width = 8, height = 6) {
  ggsave(file.path(FIGURE_DIR, filename), plot = plot, width = width, height = height)
}

assert_expected <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}
