# Programmatic Seurat/DESeq2 overlap and functional enrichment for LUAD rMos.

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

seurat_res <- read.csv(file.path(TABLE_DIR, "LUAD_rMos_FindMarkers_all.csv"), check.names = FALSE)
deseq_res <- read.csv(file.path(TABLE_DIR, "LUAD_rMos_DESeq2_significant.csv"), check.names = FALSE)

seurat_sig <- seurat_res %>% filter(!is.na(p_val_adj), p_val_adj < 0.05)
deseq_sig <- deseq_res %>% filter(!is.na(padj), padj < 0.05)

overlap_genes <- intersect(seurat_sig$gene, deseq_sig$gene)
union_genes <- union(seurat_sig$gene, deseq_sig$gene)

overlap_tbl <- seurat_sig %>%
  filter(gene %in% overlap_genes) %>%
  arrange(desc(avg_log2FC))
write.csv(overlap_tbl, file.path(TABLE_DIR, "LUAD_rMos_overlapping_DEGs.csv"), row.names = FALSE)

summary_tbl <- data.frame(
  analysis = c("Seurat", "DESeq2", "Overlap", "Union", "Overlap_percent_of_union"),
  value = c(nrow(seurat_sig), nrow(deseq_sig), length(overlap_genes), length(union_genes),
            100 * length(overlap_genes) / length(union_genes))
)
write.csv(summary_tbl, file.path(TABLE_DIR, "LUAD_rMos_DEG_overlap_summary.csv"), row.names = FALSE)

# Directionality for downstream analysis is taken from Seurat avg_log2FC.
enrichment_input <- overlap_tbl %>%
  filter(abs(avg_log2FC) > 0.5) %>%
  mutate(direction = if_else(avg_log2FC > 0, "up", "down"))
write.csv(enrichment_input, file.path(TABLE_DIR, "LUAD_rMos_overlap_absFC_gt0.5.csv"), row.names = FALSE)

id_map <- bitr(unique(enrichment_input$gene),
               fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
enrich_df <- enrichment_input %>%
  inner_join(id_map, by = c("gene" = "SYMBOL")) %>%
  distinct(gene, ENTREZID, direction, .keep_all = TRUE)

gene_clusters <- split(enrich_df$ENTREZID, enrich_df$direction)

kegg <- compareCluster(
  geneCluster = gene_clusters,
  fun = "enrichKEGG",
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
kegg <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
kegg_tbl <- as.data.frame(kegg)
write.csv(kegg_tbl, file.path(TABLE_DIR, "LUAD_rMos_KEGG_enrichment.csv"), row.names = FALSE)

run_go <- function(ontology) {
  x <- compareCluster(
    geneCluster = gene_clusters,
    fun = "enrichGO",
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = ontology,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  as.data.frame(x)
}

go_bp <- run_go("BP")
go_cc <- run_go("CC")
go_mf <- run_go("MF")
write.csv(go_bp, file.path(TABLE_DIR, "LUAD_rMos_GO_BP_enrichment.csv"), row.names = FALSE)
write.csv(go_cc, file.path(TABLE_DIR, "LUAD_rMos_GO_CC_enrichment.csv"), row.names = FALSE)
write.csv(go_mf, file.path(TABLE_DIR, "LUAD_rMos_GO_MF_enrichment.csv"), row.names = FALSE)

# Top five terms per direction by q-value for manuscript panels.
top5 <- function(x) {
  x %>% filter(!is.na(qvalue), qvalue < 0.05) %>%
    group_by(Cluster) %>%
    slice_min(order_by = qvalue, n = 5, with_ties = FALSE) %>%
    ungroup()
}
write.csv(top5(kegg_tbl), file.path(TABLE_DIR, "LUAD_rMos_KEGG_top5.csv"), row.names = FALSE)
write.csv(top5(go_bp), file.path(TABLE_DIR, "LUAD_rMos_GO_BP_top5.csv"), row.names = FALSE)
write.csv(top5(go_cc), file.path(TABLE_DIR, "LUAD_rMos_GO_CC_top5.csv"), row.names = FALSE)
write.csv(top5(go_mf), file.path(TABLE_DIR, "LUAD_rMos_GO_MF_top5.csv"), row.names = FALSE)

# Heatmap of overlapping genes using DESeq2 variance-stabilized pseudobulk data.
vst_mat <- readRDS(file.path(OBJECT_DIR, "LUAD_rMos_DESeq2_vst_matrix.rds"))
heat_mat <- vst_mat[intersect(overlap_genes, rownames(vst_mat)), , drop = FALSE]
heat_mat <- t(scale(t(heat_mat)))
heat_mat[!is.finite(heat_mat)] <- 0

pdf(file.path(FIGURE_DIR, "LUAD_rMos_overlapping_DEGs_heatmap.pdf"), width = 9, height = 12)
draw(Heatmap(
  heat_mat,
  name = "Z-score",
  show_row_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  column_title = "Overlapping Seurat/DESeq2 rMos DEGs"
))
dev.off()

cat("Seurat significant DEGs:", nrow(seurat_sig), "\n")
cat("DESeq2 significant DEGs:", nrow(deseq_sig), "\n")
cat("Overlapping genes:", length(overlap_genes), "\n")
cat("Filtered overlapping genes (|Seurat avg_log2FC| > 0.5):", nrow(enrichment_input), "\n")
cat("Up:", sum(enrichment_input$direction == "up"), " Down:", sum(enrichment_input$direction == "down"), "\n")
