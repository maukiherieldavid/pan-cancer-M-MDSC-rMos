# Melanoma independent signature evaluation and rMo-like analysis
# Dataset: GSE123139

source(file.path("scripts", "00_setup.R"))

mel_counts <- Read10X(MELANOMA_10X_DIR, gene.column = 1)
mel <- CreateSeuratObject(
  counts = mel_counts,
  project = "melanoma",
  min.cells = 3,
  min.features = 200
)
rm(mel_counts)

# QC metrics were inspected; no additional mitochondrial/count threshold was
# applied to the final melanoma object.
mel <- add_qc_metrics(mel)
mel <- normalize_pca(mel)
mel <- FindNeighbors(mel, dims = 1:20)
mel <- FindClusters(mel, resolution = 0.4)
mel <- RunUMAP(mel, dims = 1:20)

mel_metadata <- read.csv(MELANOMA_METADATA, check.names = FALSE)
meta <- mel@meta.data %>%
  rownames_to_column("cell.id") %>%
  left_join(mel_metadata, by = "cell.id") %>%
  column_to_rownames("cell.id")
mel <- AddMetaData(mel, metadata = meta)

assert_expected(ncol(mel) == 47681,
                "Melanoma cell count differs from the 47,681-cell dataset used in the manuscript.")

# Main cell-group and 15-population annotations.
Idents(mel) <- "RNA_snn_res.0.4"
cell_group_map <- setNames(
  c("Lymphoid(imm)", "Lymphoid(imm)", "Lymphoid(imm)", "Lymphoid(imm)", "Myeloid(imm)",
    "Lymphoid(imm)", "Myeloid(imm)", "Lymphoid(imm)", "Lymphoid(imm)", "Lymphoid(imm)",
    "Myeloid(imm)", "Tumor", "Myeloid(imm)", "Myeloid(imm)", "Myeloid(imm)"),
  as.character(0:14)
)
mel <- rename_idents(mel, cell_group_map)
mel$cell_group <- Idents(mel)

Idents(mel) <- "RNA_snn_res.0.4"
cell_type_map <- setNames(
  c("CD8 Tcells", "Immature NK", "B cells", "NK cells", "Monocytes", "Tcells cl", "TAMs",
    "Proliferating Tcells", "Tregs", "Memory Bcells", "DCs", "Tumor cells", "pDCs",
    "Undetermined", "Macro cl"),
  as.character(0:14)
)
mel <- rename_idents(mel, cell_type_map)
mel$cellType <- Idents(mel)

# Myeloid compartment and the 18-cluster representation used for signature evaluation.
Idents(mel) <- "cell_group"
myeloid <- subset(mel, idents = "Myeloid(imm)")
myeloid <- ScaleData(myeloid)
myeloid <- RunPCA(myeloid)
myeloid <- FindNeighbors(myeloid, dims = 1:20)
myeloid <- FindClusters(myeloid, resolution = 1.0)
myeloid <- RunUMAP(myeloid, dims = 1:20)
myeloid$myeloid_cluster <- Idents(myeloid)

assert_expected(ncol(myeloid) == 9637,
                "Melanoma myeloid cell count differs from the 9,637-cell compartment used in the manuscript.")

marker_panel <- unique(c(monocytic_context_genes, signature_genes, mdsc_context_genes))
p_signature <- DotPlot(myeloid, features = intersect(marker_panel, rownames(myeloid))) + RotatedAxis()
save_plot(p_signature, "melanoma_myeloid_signature_dotplot.pdf", width = 12, height = 6)

# Candidate clusters used for higher-resolution monocytic analysis. Clusters 2,
# 7 and 12 show the representative APOBEC3A/CD300E pattern reported in the main
# figure; cluster 4 is included in the broader monocytic compartment used for the
# high-resolution analysis.
Idents(myeloid) <- "myeloid_cluster"
candidate_ids <- intersect(c("2", "4", "7", "12"), levels(myeloid))
mono <- subset(myeloid, idents = candidate_ids)
mono <- ScaleData(mono)
mono <- RunPCA(mono)
mono <- FindNeighbors(mono, dims = 1:20)
mono <- FindClusters(mono, resolution = 0.2)
mono <- FindClusters(mono, resolution = 0.3)
mono <- RunUMAP(mono, dims = 1:20)

# Resolution 0.2 identifies the broader REL/CCR2 rMo-like population; resolution
# 0.3 resolves additional transcriptional heterogeneity including NFKB1 expression.
Idents(mono) <- "RNA_snn_res.0.2"
p_rel_ccr2 <- VlnPlot(mono, features = c("REL", "CCR2"), pt.size = 0, ncol = 2)
save_plot(p_rel_ccr2, "melanoma_rMos_REL_CCR2_res0.2.pdf", width = 9, height = 5)

Idents(mono) <- "RNA_snn_res.0.3"
p_nfkb1 <- VlnPlot(mono, features = c("IL1B", "NFKB1"), pt.size = 0, ncol = 2)
save_plot(p_nfkb1, "melanoma_rMos_IL1B_NFKB1_res0.3.pdf", width = 9, height = 5)

saveRDS(mel, file.path(OBJECT_DIR, "melanoma_all_cells.rds"))
saveRDS(myeloid, file.path(OBJECT_DIR, "melanoma_myeloid_18clusters.rds"))
saveRDS(mono, file.path(OBJECT_DIR, "melanoma_high_resolution_monocytes.rds"))

cat("Melanoma analysis complete.\n")
