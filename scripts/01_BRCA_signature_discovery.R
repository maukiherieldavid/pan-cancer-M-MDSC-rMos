# BRCA signature discovery (SCP1106)
# Reconstructs the analysis leading to the proposed five-gene M-MDSC signature.

source(file.path("scripts", "00_setup.R"))

# 1. Load counts and metadata ---------------------------------------------------
brca_counts <- read.csv(BRCA_COUNTS, row.names = 1, check.names = FALSE)
brca <- CreateSeuratObject(
  counts = brca_counts,
  project = "BRCA",
  min.cells = 3,
  min.features = 200
)
rm(brca_counts)

brca <- add_qc_metrics(brca)
brca <- normalize_pca(brca)

# 2. RPCA integration across patient/sample identities -------------------------
brca[["RNA"]] <- split(brca[["RNA"]], f = brca$orig.ident)
brca <- IntegrateLayers(
  object = brca,
  method = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.rpca",
  verbose = FALSE
)
brca <- FindNeighbors(brca, reduction = "integrated.rpca", dims = 1:20)
brca <- FindClusters(brca, resolution = 0.3)
brca <- RunUMAP(
  brca,
  reduction = "integrated.rpca",
  dims = 1:20,
  reduction.name = "umap.rpca"
)
brca <- JoinLayers(brca)

brca_metadata <- read.csv(BRCA_METADATA, check.names = FALSE)
meta <- brca@meta.data %>%
  rownames_to_column("cell_id") %>%
  left_join(brca_metadata, by = "cell_id") %>%
  column_to_rownames("cell_id")
brca <- AddMetaData(brca, metadata = meta)

# Main cellular compartments at resolution 0.3.
Idents(brca) <- "RNA_snn_res.0.3"
main_map <- setNames(
  c("immune", "immune", "immune", "epithelial", "immune", "epithelial",
    "immune", "stromal", "immune", "immune", "epithelial", "immune",
    "stromal", "stromal", "stromal", "epithelial", "epithelial"),
  as.character(0:16)
)
brca <- rename_idents(brca, main_map)
brca$cell_class <- Idents(brca)

# 3. Immune-cell subclustering -------------------------------------------------
Idents(brca) <- "cell_class"
immune <- subset(brca, idents = "immune")
immune <- ScaleData(immune)
immune <- RunPCA(immune)
immune <- FindNeighbors(immune, dims = 1:20)
immune <- FindClusters(immune, resolution = 0.3)
immune <- RunUMAP(immune, dims = 1:20)

# Exclude clusters with conflicting/undetermined lineage-marker profiles.
Idents(immune) <- "RNA_snn_res.0.3"
immune <- subset(
  immune,
  idents = c("0", "1", "2", "3", "4", "5", "7", "9", "10", "11", "12", "13", "14", "16")
)

# Final immune-cell reclustering and annotation.
immune <- ScaleData(immune)
immune <- RunPCA(immune)
immune <- FindNeighbors(immune, dims = 1:20)
immune <- FindClusters(immune, resolution = 0.4)
immune <- RunUMAP(immune, dims = 1:20)
Idents(immune) <- "RNA_snn_res.0.4"

immune_map <- setNames(
  c("Macrophages", "CD8 Tcells1", "T cells1", "CD8 Tcells2", "Plasma", "DCs",
    "T cells2", "Tregs", "B cells1", "Monocytes 1", "NK cells", "Monocytes 2",
    "B cells2", "Proliferating Tcells", "Plasma cells2", "Plasma cells3", "Neutrophils"),
  as.character(0:16)
)
immune <- rename_idents(immune, immune_map)
immune$imm_cellType <- Idents(immune)

# 4. Differential-expression analyses used for candidate discovery ------------
Idents(immune) <- "imm_cellType"
all_markers <- FindAllMarkers(
  immune,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
write.csv(all_markers, file.path(TABLE_DIR, "BRCA_FindAllMarkers.csv"), row.names = FALSE)

mono_vs_macro <- FindMarkers(
  immune,
  ident.1 = c("Monocytes 1", "Monocytes 2"),
  ident.2 = "Macrophages",
  only.pos = TRUE
)
mono_vs_all <- FindMarkers(
  immune,
  ident.1 = c("Monocytes 1", "Monocytes 2"),
  ident.2 = NULL,
  only.pos = TRUE
)

mono_vs_macro_top25 <- rank_positive_markers(mono_vs_macro, n = 25)
mono_vs_all_top25 <- rank_positive_markers(mono_vs_all, n = 25)

write.csv(mono_vs_macro %>% rownames_to_column("gene"),
          file.path(TABLE_DIR, "BRCA_monocytes_vs_macrophages_all_positive.csv"), row.names = FALSE)
write.csv(mono_vs_all %>% rownames_to_column("gene"),
          file.path(TABLE_DIR, "BRCA_monocytes_vs_all_immune_all_positive.csv"), row.names = FALSE)
save_table(mono_vs_macro_top25, "BRCA_monocytes_vs_macrophages_top25.csv")
save_table(mono_vs_all_top25, "BRCA_monocytes_vs_all_immune_top25.csv")

# Candidate expression was evaluated across all annotated immune populations.
# The final five genes were retained because they showed preferential expression
# in both monocytic clusters with comparatively limited expression elsewhere.
p_signature <- VlnPlot(
  immune,
  features = c("CD14", signature_genes),
  pt.size = 0,
  ncol = 3
)
save_plot(p_signature, "BRCA_five_gene_signature_violin.pdf", width = 12, height = 8)

p_allimmune <- ggplot(mono_vs_all_top25,
  aes(x = avg_log2FC, y = reorder(gene, avg_log2FC))) +
  geom_col() +
  labs(x = "Average log2 fold change", y = NULL,
       title = "BRCA monocytic clusters vs all other immune populations") +
  theme_classic()
save_plot(p_allimmune, "BRCA_monocytes_vs_all_immune_top25.pdf", width = 8, height = 7)

saveRDS(brca, file.path(OBJECT_DIR, "BRCA_all_cells.rds"))
saveRDS(immune, file.path(OBJECT_DIR, "BRCA_immune_annotated.rds"))

cat("BRCA analysis complete.\n")
cat("Immune cells:", ncol(immune), "\n")
cat("Final signature:", paste(signature_genes, collapse = ", "), "\n")
