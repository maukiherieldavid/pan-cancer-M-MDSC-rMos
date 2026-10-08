# Recreate core computational panels from saved analysis objects/tables.
# Microscopy panels and final multi-panel layout assembly are not generated here.

source(file.path("scripts", "00_setup.R"))

# Helper -----------------------------------------------------------------------
plot_if_reduction <- function(object, reduction, group.by, filename, width = 7, height = 6) {
  if (reduction %in% Reductions(object) && group.by %in% colnames(object@meta.data)) {
    p <- DimPlot(object, reduction = reduction, group.by = group.by, label = TRUE, repel = TRUE)
    save_plot(p, filename, width, height)
  }
}

# Figure 1: dataset-level cell compartments ------------------------------------
brca <- readRDS(file.path(OBJECT_DIR, "BRCA_all_cells.rds"))
crc <- readRDS(file.path(OBJECT_DIR, "CRC_all_cells.rds"))
mel <- readRDS(file.path(OBJECT_DIR, "melanoma_all_cells.rds"))
luad <- readRDS(file.path(OBJECT_DIR, "LUAD_all_cells_harmony.rds"))

plot_if_reduction(brca, "umap.rpca", "cell_class", "Fig1_BRCA_compartments_UMAP.pdf")
plot_if_reduction(crc, "umap", "cell_group", "Fig1_CRC_compartments_UMAP.pdf")
plot_if_reduction(mel, "umap", "cell_group", "Fig1_melanoma_compartments_UMAP.pdf")
plot_if_reduction(luad, "umap.harmony", "cell_group", "Fig1_LUAD_compartments_UMAP.pdf")

# Figure 2: BRCA immune annotation and signature -------------------------------
brca_immune <- readRDS(file.path(OBJECT_DIR, "BRCA_immune_annotated.rds"))
plot_if_reduction(brca_immune, "umap", "imm_cellType", "Fig2a_BRCA_immune_UMAP.pdf", 9, 7)
Idents(brca_immune) <- "imm_cellType"
p_br <- VlnPlot(brca_immune, features = c("CD14", signature_genes), pt.size = 0, ncol = 3)
save_plot(p_br, "Fig2e_BRCA_signature_violin.pdf", width = 12, height = 8)

# Figure 4 / Figure S2: independent cross-cancer evaluation --------------------
crc_immune <- readRDS(file.path(OBJECT_DIR, "CRC_immune_annotated.rds"))
crc_myeloid <- readRDS(file.path(OBJECT_DIR, "CRC_myeloid_signature_clusters.rds"))
mel_myeloid <- readRDS(file.path(OBJECT_DIR, "melanoma_myeloid_18clusters.rds"))
luad_immune <- readRDS(file.path(OBJECT_DIR, "LUAD_immune_annotated.rds"))
luad_myeloid <- readRDS(file.path(OBJECT_DIR, "LUAD_myeloid_clusters.rds"))

plot_if_reduction(crc_immune, "umap", "imm_cellType", "Fig4a_CRC_immune_UMAP.pdf")
plot_if_reduction(crc_myeloid, "umap", "myeloid_cluster", "Fig4b_CRC_myeloid_UMAP.pdf")
plot_if_reduction(mel, "umap", "cellType", "Fig4c_melanoma_annotated_UMAP.pdf")
plot_if_reduction(mel_myeloid, "umap", "myeloid_cluster", "Fig4d_melanoma_myeloid_UMAP.pdf")
plot_if_reduction(luad_immune, "umap.harmony", "cellType", "Fig4e_LUAD_immune_UMAP.pdf")
plot_if_reduction(luad_myeloid, "umap.harmony", "myeloid_cluster", "Fig4f_LUAD_myeloid_UMAP.pdf")

for (spec in list(
  list(name = "CRC", object = crc_myeloid, group = "myeloid_cluster"),
  list(name = "melanoma", object = mel_myeloid, group = "myeloid_cluster"),
  list(name = "LUAD", object = luad_myeloid, group = "myeloid_cluster")
)) {
  obj <- spec$object
  panel <- intersect(unique(c("CD14", "FCN1", "TGFB1", signature_genes)), rownames(obj))
  p <- DotPlot(obj, features = panel, group.by = spec$group) + RotatedAxis()
  save_plot(p, paste0("Fig4_FigS2_", spec$name, "_candidate_markers.pdf"), width = 10, height = 6)
}

# Figure 5: high-resolution monocytic populations ------------------------------
crc_mono <- readRDS(file.path(OBJECT_DIR, "CRC_high_resolution_monocytes.rds"))
mel_mono <- readRDS(file.path(OBJECT_DIR, "melanoma_high_resolution_monocytes.rds"))
luad_mono <- readRDS(file.path(OBJECT_DIR, "LUAD_high_resolution_monocytes.rds"))

plot_if_reduction(crc_mono, "umap", "highres_cluster", "Fig5a_CRC_highres_UMAP.pdf")
plot_if_reduction(mel_mono, "umap", "RNA_snn_res.0.2", "Fig5c_melanoma_res0.2_UMAP.pdf")
plot_if_reduction(mel_mono, "umap", "RNA_snn_res.0.3", "Fig5c_melanoma_res0.3_UMAP.pdf")
plot_if_reduction(luad_mono, "umap.harmony", "mono_subsets", "Fig5b_Fig6a_LUAD_highres_UMAP.pdf")

for (spec in list(
  list(name = "CRC", object = crc_mono, group = "highres_cluster"),
  list(name = "melanoma_res0.2", object = mel_mono, group = "RNA_snn_res.0.2"),
  list(name = "LUAD", object = luad_mono, group = "mono_subsets")
)) {
  p <- DotPlot(spec$object,
               features = intersect(c("REL", "CCR2", "IL1B", "NFKB1"), rownames(spec$object)),
               group.by = spec$group) + RotatedAxis()
  save_plot(p, paste0("Fig5_", spec$name, "_rMo_markers.pdf"), width = 9, height = 5)
}

# Figure 6b: LUAD monocytic subset markers -------------------------------------
p_markers <- DotPlot(
  luad_mono,
  features = intersect(c("CXCL8", "HMOX1", "SPP1", "RETN", "CXCL9", "REL", "CCR2", "NFKB1"), rownames(luad_mono)),
  group.by = "mono_subsets"
) + RotatedAxis()
save_plot(p_markers, "Fig6b_LUAD_monocytic_subset_markers.pdf", width = 9, height = 5)

cat("Computational figure-panel export complete.\n")
