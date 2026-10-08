# Programmatic Seurat/DESeq2 overlap and functional enrichment for LUAD rMos.
# The enrichment workflow follows the analysis used for Tables S2-S5 and
# Figure 6f-g / Figure S5.

source(file.path("scripts", "00_setup.R"))

for (pkg in c("clusterProfiler", "org.Hs.eg.db", "ComplexHeatmap", "circlize")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Install required package: ", pkg)
}

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ComplexHeatmap)
  library(circlize)
})

# -----------------------------------------------------------------------------
# 1. Identify the overlapping Seurat and DESeq2 DEG sets
# -----------------------------------------------------------------------------

seurat_res <- read.csv(
  file.path(TABLE_DIR, "LUAD_rMos_FindMarkers_all.csv"),
  check.names = FALSE
)

deseq_res <- read.csv(
  file.path(TABLE_DIR, "LUAD_rMos_DESeq2_significant.csv"),
  check.names = FALSE
)

seurat_sig <- seurat_res %>%
  filter(!is.na(p_val_adj), p_val_adj < 0.05)

deseq_sig <- deseq_res %>%
  filter(!is.na(padj), padj < 0.05)

overlap_genes <- intersect(seurat_sig$gene, deseq_sig$gene)
union_genes <- union(seurat_sig$gene, deseq_sig$gene)

overlap_tbl <- seurat_sig %>%
  filter(gene %in% overlap_genes) %>%
  arrange(desc(avg_log2FC))

write.csv(
  overlap_tbl,
  file.path(TABLE_DIR, "LUAD_rMos_overlapping_DEGs.csv"),
  row.names = FALSE
)

summary_tbl <- data.frame(
  analysis = c("Seurat", "DESeq2", "Overlap", "Union", "Overlap_percent_of_union"),
  value = c(
    nrow(seurat_sig),
    nrow(deseq_sig),
    length(overlap_genes),
    length(union_genes),
    100 * length(overlap_genes) / length(union_genes)
  )
)

write.csv(
  summary_tbl,
  file.path(TABLE_DIR, "LUAD_rMos_DEG_overlap_summary.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 2. Prepare the shared DEG set used for enrichment
# -----------------------------------------------------------------------------
# Directionality follows the Seurat avg_log2FC values reported in the study.

enrichment_input <- overlap_tbl %>%
  filter(abs(avg_log2FC) > 0.5) %>%
  mutate(group = if_else(avg_log2FC > 0.5, "up", "down"))

write.csv(
  enrichment_input,
  file.path(TABLE_DIR, "LUAD_rMos_overlap_absFC_gt0.5.csv"),
  row.names = FALSE
)

Gene_ID <- bitr(
  enrichment_input$gene,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

enrichment_data <- merge(
  Gene_ID,
  enrichment_input[, c("gene", "group")],
  by.x = "SYMBOL",
  by.y = "gene"
) %>%
  distinct(SYMBOL, ENTREZID, group)

options(timeout = 3600)

# Utility: standardize the grouping column returned by compareCluster.
standardize_group <- function(x) {
  if (!"group" %in% colnames(x) && "Cluster" %in% colnames(x)) {
    x$group <- x$Cluster
  }
  x
}

# Utility: retain the five most significant terms per direction by q-value.
top5_by_qvalue <- function(x) {
  x <- standardize_group(x)
  x %>%
    filter(!is.na(qvalue), qvalue < 0.05) %>%
    group_by(group) %>%
    slice_min(order_by = qvalue, n = 5, with_ties = FALSE) %>%
    ungroup()
}

# -----------------------------------------------------------------------------
# 3. KEGG pathway enrichment (Table S2; Figure 6f)
# -----------------------------------------------------------------------------

kegg_obj <- compareCluster(
  ENTREZID ~ group,
  data = enrichment_data,
  fun = "enrichKEGG",
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

kegg_obj <- setReadable(
  kegg_obj,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
)

kegg_tbl <- standardize_group(kegg_obj@compareClusterResult)

# Complete KEGG enrichment output used for Table S2.
write.csv(
  kegg_tbl,
  file.path(TABLE_DIR, "Table_S2_KEGG_enrichment.csv"),
  row.names = FALSE
)

kegg_top5 <- top5_by_qvalue(kegg_tbl)
write.csv(
  kegg_top5,
  file.path(TABLE_DIR, "LUAD_rMos_KEGG_top5.csv"),
  row.names = FALSE
)

# Display up to five genes per pathway, matching the manuscript visualization.
kegg_plot <- kegg_top5 %>%
  mutate(
    group = factor(group, levels = c("up", "down")),
    gene_label = vapply(
      strsplit(geneID, "/", fixed = TRUE),
      function(x) paste(head(x, 5), collapse = "/"),
      character(1)
    )
  ) %>%
  arrange(group, qvalue) %>%
  mutate(Description = factor(Description, levels = rev(unique(Description))))

p_kegg <- ggplot(
  kegg_plot,
  aes(x = -log10(qvalue), y = Description, fill = group)
) +
  geom_col(width = 0.5) +
  geom_text(
    aes(x = 0.1, label = Description),
    size = 3.5,
    hjust = 0
  ) +
  geom_text(
    aes(x = 0.1, label = gene_label, colour = group),
    size = 4,
    fontface = "italic",
    hjust = 0,
    vjust = 2.3
  ) +
  scale_fill_manual(values = c("up" = "#CB5640", "down" = "#65B0C6")) +
  scale_colour_manual(values = c("up" = "#CB5640", "down" = "#65B0C6")) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(
    title = "KEGG Pathway Enrichment",
    x = "-log10(q-value)",
    y = "Pathways"
  ) +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_text(colour = "black", size = 12),
    axis.line = element_line(colour = "black", linewidth = 0.5),
    axis.text.x = element_text(colour = "black", size = 10),
    axis.ticks.x = element_line(colour = "black"),
    axis.title.x = element_text(colour = "black", size = 12),
    legend.position = "none"
  )

save_plot(p_kegg, "Fig6f_KEGG_enrichment.pdf", width = 5, height = 5)
ggsave(
  file.path(FIGURE_DIR, "Fig6f_KEGG_enrichment.png"),
  plot = p_kegg,
  width = 5,
  height = 5,
  dpi = 300
)

# -----------------------------------------------------------------------------
# 4. Gene Ontology enrichment (Tables S3-S5; Figure 6g and Figure S5)
# -----------------------------------------------------------------------------

run_go <- function(ontology) {
  go_obj <- compareCluster(
    ENTREZID ~ group,
    data = enrichment_data,
    fun = "enrichGO",
    OrgDb = org.Hs.eg.db,
    ont = ontology,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  standardize_group(go_obj@compareClusterResult)
}

go_bp <- run_go("BP")
go_cc <- run_go("CC")
go_mf <- run_go("MF")

write.csv(go_bp, file.path(TABLE_DIR, "Table_S3_GO_BP_enrichment.csv"), row.names = FALSE)
write.csv(go_cc, file.path(TABLE_DIR, "Table_S4_GO_CC_enrichment.csv"), row.names = FALSE)
write.csv(go_mf, file.path(TABLE_DIR, "Table_S5_GO_MF_enrichment.csv"), row.names = FALSE)

plot_go_split <- function(go_tbl, ontology_title, file_stub) {
  go_top <- top5_by_qvalue(go_tbl) %>%
    mutate(
      log10_qvalue = if_else(group == "up", -log10(qvalue), log10(qvalue)),
      short_description = if_else(
        nchar(Description) > 35,
        paste0(substr(Description, 1, 32), "..."),
        Description
      )
    )

  go_down <- go_top %>% filter(group == "down") %>% arrange(qvalue)
  go_up <- go_top %>% filter(group == "up") %>% arrange(qvalue)
  go_plot <- bind_rows(go_down, go_up)

  go_plot$Description <- factor(
    go_plot$Description,
    levels = go_plot$Description[order(go_plot$log10_qvalue)]
  )

  x_max <- max(abs(go_plot$log10_qvalue), na.rm = TRUE) + 0.5

  p <- ggplot(
    go_plot,
    aes(x = log10_qvalue, y = Description, fill = group)
  ) +
    geom_col(width = 0.6) +
    scale_x_continuous(limits = c(-x_max, x_max)) +
    scale_fill_manual(values = c("up" = "#CB5640", "down" = "#65B0C6")) +
    geom_text(
      data = go_plot %>% filter(group == "up"),
      aes(x = -0.1, label = short_description, colour = group),
      size = 4.5,
      hjust = 1,
      fontface = "bold"
    ) +
    geom_text(
      data = go_plot %>% filter(group == "down"),
      aes(x = 0.1, label = short_description, colour = group),
      size = 4.5,
      hjust = 0,
      fontface = "bold"
    ) +
    scale_colour_manual(values = c("up" = "#CB5640", "down" = "#65B0C6")) +
    labs(
      x = expression("log"[10] * "(q-value)"),
      y = NULL,
      title = ontology_title
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 12, face = "bold"),
      axis.text.y = element_blank(),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "none",
      plot.margin = margin(t = 20, r = 40, b = 20, l = 40)
    ) +
    geom_vline(xintercept = 0, linetype = "solid", colour = "black", alpha = 0.3)

  save_plot(p, paste0(file_stub, ".pdf"), width = 10, height = 5)
  ggsave(
    file.path(FIGURE_DIR, paste0(file_stub, ".png")),
    plot = p,
    width = 10,
    height = 5,
    dpi = 300
  )

  invisible(go_top)
}

bp_top5 <- plot_go_split(
  go_bp,
  "Gene Ontology Biological Process Enrichment",
  "Fig6g_GO_BP_enrichment"
)

cc_top5 <- plot_go_split(
  go_cc,
  "Gene Ontology Cellular Component Enrichment",
  "FigS5a_GO_CC_enrichment"
)

mf_top5 <- plot_go_split(
  go_mf,
  "Gene Ontology Molecular Function Enrichment",
  "FigS5b_GO_MF_enrichment"
)

write.csv(bp_top5, file.path(TABLE_DIR, "LUAD_rMos_GO_BP_top5.csv"), row.names = FALSE)
write.csv(cc_top5, file.path(TABLE_DIR, "LUAD_rMos_GO_CC_top5.csv"), row.names = FALSE)
write.csv(mf_top5, file.path(TABLE_DIR, "LUAD_rMos_GO_MF_top5.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 5. Hierarchical clustering of overlapping DEGs (Figure 6e)
# -----------------------------------------------------------------------------

vst_mat <- readRDS(file.path(OBJECT_DIR, "LUAD_rMos_DESeq2_vst_matrix.rds"))
heat_mat <- vst_mat[intersect(overlap_genes, rownames(vst_mat)), , drop = FALSE]
heat_mat <- t(scale(t(heat_mat)))
heat_mat[!is.finite(heat_mat)] <- 0

pdf(
  file.path(FIGURE_DIR, "Fig6e_LUAD_rMos_overlapping_DEGs_heatmap.pdf"),
  width = 9,
  height = 12
)
draw(
  Heatmap(
    heat_mat,
    name = "Z-score",
    show_row_names = FALSE,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    column_title = "Overlapping Seurat/DESeq2 rMos DEGs"
  )
)
dev.off()

# -----------------------------------------------------------------------------
# 6. Analysis summary
# -----------------------------------------------------------------------------

cat("Seurat significant DEGs:", nrow(seurat_sig), "\n")
cat("DESeq2 significant DEGs:", nrow(deseq_sig), "\n")
cat("Overlapping genes:", length(overlap_genes), "\n")
cat("Union genes:", length(union_genes), "\n")
cat("Overlap (% of union):", round(100 * length(overlap_genes) / length(union_genes), 1), "\n")
cat("Filtered overlapping genes (|Seurat avg_log2FC| > 0.5):", nrow(enrichment_input), "\n")
cat("Up:", sum(enrichment_input$group == "up"), "Down:", sum(enrichment_input$group == "down"), "\n")
cat("KEGG enriched pathways:", nrow(kegg_tbl), "\n")
cat("GO BP enriched terms:", nrow(go_bp), "\n")
cat("GO CC enriched terms:", nrow(go_cc), "\n")
cat("GO MF enriched terms:", nrow(go_mf), "\n")
