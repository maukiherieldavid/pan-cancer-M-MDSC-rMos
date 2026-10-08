# LUAD preprocessing, independent signature evaluation, and rMos identification
# Dataset: GSE131907

source(file.path("scripts", "00_setup.R"))

normal_samples <- c("p018n", "p019n", "p027n", "p028n", "p029n", "p030n", "p031n", "p032n", "p033n", "p034n")
tumour_samples <- c("p018t", "p019t", "p023t", "p024t", "p027t", "p030t", "p031t", "p032t", "p033t", "p034t")
all_samples <- c(normal_samples, tumour_samples)

sample_dir <- function(sample_id) {
  parent <- if (str_ends(sample_id, "n")) LUAD_NORMAL_DIR else LUAD_TUMOR_DIR
  file.path(parent, paste0("filtered_feature_bc_matrix.", sample_id))
}

# 1. Read each 10x sample and merge -------------------------------------------
luad_list <- lapply(all_samples, function(sid) {
  m <- Read10X(sample_dir(sid))
  CreateSeuratObject(counts = m, min.cells = 3, min.features = 200, project = "LUAD")
})
names(luad_list) <- all_samples

luad <- merge(
  x = luad_list[[1]],
  y = luad_list[-1],
  add.cell.ids = all_samples,
  project = "LUAD"
)
luad <- JoinLayers(luad)
rm(luad_list)

luad <- add_qc_metrics(luad)
luad <- subset(luad, subset = nCount_RNA < 60000)

# Derive sample/patient/tissue metadata directly from the prefixed cell names.
cell_ids <- colnames(luad)
sample_id <- str_extract(cell_ids, "^p\\d{3}[nt]")
luad$sample_id <- sample_id
luad$patient_id <- str_remove(sample_id, "[nt]$")
luad$tissue_type <- ifelse(str_ends(sample_id, "t"), "Tumor", "Normal")
luad$tissue_type <- factor(luad$tissue_type, levels = c("Normal", "Tumor"))

luad <- normalize_pca(luad)

# 2. Harmony integration split by tissue type --------------------------------
luad[["RNA"]] <- split(luad[["RNA"]], f = luad$tissue_type)
luad <- IntegrateLayers(
  object = luad,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony",
  verbose = FALSE
)
luad <- JoinLayers(luad)
luad <- FindNeighbors(luad, reduction = "harmony", dims = 1:20)
luad <- FindClusters(luad, resolution = 0.1)
luad <- RunUMAP(luad, reduction = "harmony", dims = 1:20, reduction.name = "umap.harmony")

# 3. Main compartments ---------------------------------------------------------
Idents(luad) <- "RNA_snn_res.0.1"
main_map <- setNames(
  c("immune", "immune", "immune", "stromal", "epithelial", "epithelial", "immune",
    "epithelial", "epithelial", "epithelial", "epithelial", "stromal", "unknown", "immune"),
  as.character(0:13)
)
luad <- rename_idents(luad, main_map)
luad$cell_group <- Idents(luad)
luad <- subset(luad, idents = c("immune", "epithelial", "stromal"))

# 4. Immune-cell compartment ---------------------------------------------------
# The downstream immune compartment was restricted to percent.mt < 15%.
immune <- subset(luad, subset = percent.mt < 15)
Idents(immune) <- "cell_group"
immune <- subset(immune, idents = "immune")
immune <- ScaleData(immune)
immune <- FindNeighbors(immune, reduction = "harmony", dims = 1:20)
immune <- FindClusters(immune, resolution = 0.2)
immune <- RunUMAP(immune, reduction = "harmony", dims = 1:20, reduction.name = "umap.harmony")

Idents(immune) <- "RNA_snn_res.0.2"
immune_map <- setNames(
  c("Tcells cl1", "Macro", "Mon-derived macro", "NK cells", "Alveolar macro1", "Monocytes",
    "DCs", "Tcells cl2", "B/Plasma", "Mast cells", "Alveolar macro2", "Cycling macro",
    "Proliferating Tcells"),
  as.character(0:12)
)
immune <- rename_idents(immune, immune_map)
immune$cellType <- Idents(immune)

# 5. Myeloid subclustering -----------------------------------------------------
Idents(immune) <- "cellType"
myeloid_types <- c("Macro", "Mon-derived macro", "Alveolar macro1", "Monocytes", "DCs", "Mast cells", "Alveolar macro2", "Cycling macro")
myeloid <- subset(immune, idents = myeloid_types)
myeloid <- ScaleData(myeloid)
myeloid <- FindNeighbors(myeloid, reduction = "harmony", dims = 1:20)
myeloid <- FindClusters(myeloid, resolution = 0.4)
myeloid <- RunUMAP(myeloid, reduction = "harmony", dims = 1:20, reduction.name = "umap.harmony")
myeloid$myeloid_cluster <- Idents(myeloid)

marker_panel <- unique(c(monocytic_context_genes, signature_genes, mdsc_context_genes))
p_signature <- DotPlot(
  myeloid,
  features = intersect(marker_panel, rownames(myeloid)),
  group.by = "myeloid_cluster"
) + RotatedAxis()
save_plot(p_signature, "LUAD_myeloid_signature_dotplot.pdf", width = 12, height = 6)

# Candidate monocytic compartment used for high-resolution analysis.
Idents(myeloid) <- "myeloid_cluster"
candidate_ids <- intersect(c("4", "5", "6", "7"), levels(myeloid))
mono <- subset(myeloid, idents = candidate_ids)
mono <- ScaleData(mono)
mono <- FindNeighbors(mono, reduction = "harmony", dims = 1:20)
mono <- FindClusters(mono, resolution = 0.3)
mono <- RunUMAP(mono, reduction = "harmony", dims = 1:20, reduction.name = "umap.harmony")

# Combine the two REL/CCR2-associated rMos subclusters into a single rMos state,
# retaining four additional monocytic states.
Idents(mono) <- "RNA_snn_res.0.3"
cluster_levels <- levels(mono)
assert_expected(length(cluster_levels) == 6,
                "LUAD high-resolution clustering did not produce the six source clusters expected before rMos consolidation.")
merge_map <- setNames(c("0", "1", "2", "3", "4", "1"), cluster_levels)
mono <- rename_idents(mono, merge_map)
mono$RNA_snn_res.0.3b <- Idents(mono)

Idents(mono) <- "RNA_snn_res.0.3b"
subset_map <- setNames(
  c("RETN+ monocytes", "rMos", "HMOX1+ monocytes", "CXCL9+ monocytes", "CXCL8+ monocytes"),
  levels(mono)
)
mono <- rename_idents(mono, subset_map)
mono$mono_subsets <- Idents(mono)

p_rmo <- DotPlot(mono, features = intersect(c("REL", "CCR2", "IL1B", "NFKB1", "SPP1"), rownames(mono))) + RotatedAxis()
save_plot(p_rmo, "LUAD_high_resolution_rMos_markers.pdf", width = 9, height = 5)

saveRDS(luad, file.path(OBJECT_DIR, "LUAD_all_cells_harmony.rds"))
saveRDS(immune, file.path(OBJECT_DIR, "LUAD_immune_annotated.rds"))
saveRDS(myeloid, file.path(OBJECT_DIR, "LUAD_myeloid_clusters.rds"))
saveRDS(mono, file.path(OBJECT_DIR, "LUAD_high_resolution_monocytes.rds"))

cat("LUAD signature/rMos analysis complete.\n")
