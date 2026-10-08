# CRC independent signature evaluation and high-resolution rMo-like analysis
# Dataset: GSE178341

source(file.path("scripts", "00_setup.R"))

# 1. Load the complete CRC count matrix ----------------------------------------
crc_counts <- Read10X_h5(CRC_COUNTS_H5, use.names = TRUE, unique.features = TRUE)
crc <- CreateSeuratObject(
  counts = crc_counts,
  project = "CRC",
  min.cells = 3,
  min.features = 200
)
rm(crc_counts)

# QC metrics were inspected; no additional mitochondrial/count threshold was
# applied to the final CRC atlas after initial object construction.
crc <- add_qc_metrics(crc)
crc <- normalize_pca(crc)
crc <- FindNeighbors(crc, dims = 1:20)
crc <- FindClusters(crc, resolution = 0.3)
crc <- RunUMAP(crc, dims = 1:20)

crc_metadata <- read.csv(CRC_METADATA, check.names = FALSE)
# The public metadata table is joined by cellID.
meta <- crc@meta.data %>%
  rownames_to_column("cellID") %>%
  left_join(crc_metadata, by = "cellID") %>%
  column_to_rownames("cellID")
crc <- AddMetaData(crc, metadata = meta)

assert_expected(ncol(crc) == 370114,
                "CRC cell count differs from the 370,114-cell dataset used in the manuscript.")

# 2. Main compartments ---------------------------------------------------------
Idents(crc) <- "RNA_snn_res.0.3"
main_map <- setNames(
  c("epithelial", "immune", "immune", "epithelial", "immune", "epithelial",
    "immune", "immune", "epithelial", "immune", "epithelial", "stromal",
    "epithelial", "epithelial", "stromal", "immune", "epithelial", "immune",
    "immune", "immune", "epithelial", "stromal", "immune", "epithelial"),
  as.character(0:23)
)
crc <- rename_idents(crc, main_map)
crc$cell_group <- Idents(crc)

# 3. Immune subclustering ------------------------------------------------------
Idents(crc) <- "cell_group"
immune <- subset(crc, idents = "immune")
immune <- ScaleData(immune)
immune <- RunPCA(immune)
immune <- FindNeighbors(immune, dims = 1:20)
immune <- FindClusters(immune, resolution = 0.3)
immune <- RunUMAP(immune, dims = 1:20)
Idents(immune) <- "RNA_snn_res.0.3"

immune_map <- setNames(
  c("NK/T cells", "Tregs", "Neutrophils", "mregDCs", "Immature NK",
    "Cycling Bcells", "pDCs", "Immature Bcells", "Macrophages",
    "Matured Bcell1", "Monocytes", "Plasma cells", "Matured Bcell2",
    "Mast cells", "Proliferating Tcells"),
  as.character(0:14)
)
immune <- rename_idents(immune, immune_map)
immune$imm_cellType <- Idents(immune)

# Retain a broad myeloid grouping based on the immune clusters used in the
# myeloid analysis. Clusters with conflicting lineage-marker profiles are excluded
# before final annotation and downstream analysis.
Idents(immune) <- "RNA_snn_res.0.3"
myeloid_source_ids <- c("2", "3", "6", "8", "10", "13")
myeloid <- subset(immune, idents = myeloid_source_ids)
myeloid <- ScaleData(myeloid)
myeloid <- RunPCA(myeloid)
myeloid <- FindNeighbors(myeloid, dims = 1:20)
myeloid <- FindClusters(myeloid, resolution = 0.1)
myeloid <- RunUMAP(myeloid, dims = 1:20)

# Exclude conflicting-lineage clusters and recluster.
Idents(myeloid) <- "RNA_snn_res.0.1"
myeloid <- subset(myeloid, idents = c("0", "1", "2", "3", "4", "5", "6", "8"))
myeloid <- ScaleData(myeloid)
myeloid <- RunPCA(myeloid)
myeloid <- FindNeighbors(myeloid, dims = 1:18)
myeloid <- FindClusters(myeloid, resolution = 0.2)
myeloid <- RunUMAP(myeloid, dims = 1:18)

# Reassess lineage-marker profiles, retain the myeloid clusters, and recluster.
Idents(myeloid) <- "RNA_snn_res.0.2"
myeloid <- subset(myeloid, idents = c("0", "1", "2", "3", "4", "5", "7", "8", "9"))
myeloid <- ScaleData(myeloid)
myeloid <- RunPCA(myeloid)
myeloid <- FindNeighbors(myeloid, dims = 1:30)
myeloid <- FindClusters(myeloid, resolution = 0.3)
myeloid <- RunUMAP(myeloid, dims = 1:30)

# This numbered myeloid representation is used for candidate-marker evaluation.
myeloid$myeloid_cluster <- Idents(myeloid)
myeloid_signature <- myeloid

# 4. Candidate M-MDSC-associated marker evaluation ----------------------------
Idents(myeloid_signature) <- "myeloid_cluster"
marker_panel <- unique(c(monocytic_context_genes, signature_genes, mdsc_context_genes))
p_signature <- DotPlot(myeloid_signature, features = intersect(marker_panel, rownames(myeloid_signature))) +
  RotatedAxis()
save_plot(p_signature, "CRC_myeloid_signature_dotplot.pdf", width = 11, height = 6)

# Representative candidate CRC clusters highlighted in the manuscript.
candidate_crc <- intersect(c("2", "10"), levels(myeloid_signature))
if (length(candidate_crc) > 0) {
  candidate_cells <- subset(myeloid_signature, idents = candidate_crc)
  saveRDS(candidate_cells, file.path(OBJECT_DIR, "CRC_candidate_M_MDSC_clusters.rds"))
}

# Final marker-based cleanup and broad myeloid annotation.
Idents(myeloid) <- "RNA_snn_res.0.3"
final_retained_ids <- c("0", "1", "2", "5", "6", "7", "8", "10", "11", "12")
final_retained_ids <- intersect(final_retained_ids, levels(myeloid))
myeloid_final <- subset(myeloid, idents = final_retained_ids)
myeloid_final <- ScaleData(myeloid_final)
myeloid_final <- RunPCA(myeloid_final)
myeloid_final <- FindNeighbors(myeloid_final, dims = 1:12)
myeloid_final <- FindClusters(myeloid_final, resolution = 0.2)
myeloid_final <- RunUMAP(myeloid_final, dims = 1:12)

Idents(myeloid_final) <- "RNA_snn_res.0.2"
assert_expected(length(levels(myeloid_final)) == 10,
                "CRC final myeloid clustering did not produce the ten source clusters expected before consolidation.")

# Consolidate closely related macrophage clusters and assign broad myeloid classes.
myeloid_type_map <- c(
  "0" = "Macrophages", "1" = "Monocytes", "2" = "cDCs", "3" = "Neutrophils",
  "4" = "Cycling Macro", "5" = "mregDCs", "6" = "pDCs",
  "7" = "Macrophages", "8" = "Macrophages", "9" = "Mast cells"
)
myeloid_final <- rename_idents(myeloid_final, myeloid_type_map)
myeloid_final$myeloid_cell_type <- Idents(myeloid_final)

# 5. High-resolution monocytic analysis ---------------------------------------
Idents(myeloid_final) <- "myeloid_cell_type"
mono <- subset(myeloid_final, idents = "Monocytes")
mono <- ScaleData(mono)
mono <- RunPCA(mono)
mono <- FindNeighbors(mono, dims = 1:20)
mono <- FindClusters(mono, resolution = 0.2)
mono <- RunUMAP(mono, dims = 1:20)
mono$highres_cluster <- Idents(mono)

p_rmo <- VlnPlot(
  mono,
  features = intersect(c("REL", "CCR2", "IL1B", "NFKB1"), rownames(mono)),
  pt.size = 0,
  ncol = 2
)
save_plot(p_rmo, "CRC_highres_REL_CCR2_NFKB1.pdf", width = 10, height = 7)

saveRDS(myeloid_signature, file.path(OBJECT_DIR, "CRC_myeloid_signature_clusters.rds"))
saveRDS(myeloid_final, file.path(OBJECT_DIR, "CRC_myeloid_final.rds"))
saveRDS(mono, file.path(OBJECT_DIR, "CRC_high_resolution_monocytes.rds"))

saveRDS(crc, file.path(OBJECT_DIR, "CRC_all_cells.rds"))
saveRDS(immune, file.path(OBJECT_DIR, "CRC_immune_annotated.rds"))

cat("CRC analysis complete.\n")
