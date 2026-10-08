# Sample-level pseudobulk DESeq2 analysis of LUAD rMos.
# Tissue type is the design factor, matching the analysis reported in the manuscript.

source(file.path("scripts", "00_setup.R"))

for (pkg in c("DESeq2", "apeglm")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Install Bioconductor package: ", pkg)
}
suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
})

mono <- readRDS(file.path(OBJECT_DIR, "LUAD_high_resolution_monocytes.rds"))
Idents(mono) <- "mono_subsets"
rmos <- subset(mono, idents = "rMos")

counts <- GetAssayData(rmos, assay = "RNA", layer = "counts")
meta <- rmos@meta.data %>%
  rownames_to_column("cell_id") %>%
  select(cell_id, sample_id, patient_id, tissue_type)

# Sparse sample aggregation: genes x cells multiplied by a cells x samples
# incidence matrix.
sample_factor <- factor(meta$sample_id, levels = unique(meta$sample_id))
incidence <- sparse.model.matrix(~ 0 + sample_factor)
colnames(incidence) <- levels(sample_factor)
pseudobulk_counts <- counts %*% incidence
pseudobulk_counts <- round(as.matrix(pseudobulk_counts))

sample_meta <- meta %>%
  distinct(sample_id, patient_id, tissue_type) %>%
  arrange(match(sample_id, colnames(pseudobulk_counts))) %>%
  column_to_rownames("sample_id")

stopifnot(identical(rownames(sample_meta), colnames(pseudobulk_counts)))
sample_meta$tissue_type <- relevel(factor(sample_meta$tissue_type), ref = "Normal")

dds <- DESeqDataSetFromMatrix(
  countData = pseudobulk_counts,
  colData = sample_meta,
  design = ~ tissue_type
)
dds <- dds[rowSums(counts(dds)) > 0, ]
dds <- DESeq(dds)

res <- results(dds, contrast = c("tissue_type", "Tumor", "Normal"), alpha = 0.05)
coef_name <- grep("tissue_type.*Tumor.*Normal", resultsNames(dds), value = TRUE)
if (length(coef_name) == 1) {
  res_shrunk <- lfcShrink(dds, coef = coef_name, type = "apeglm")
} else {
  res_shrunk <- res
}

res_tbl <- as.data.frame(res_shrunk) %>% rownames_to_column("gene")
sig_res <- res_tbl %>% filter(!is.na(padj), padj < 0.05)

write.csv(res_tbl, file.path(TABLE_DIR, "LUAD_rMos_DESeq2_all.csv"), row.names = FALSE)
write.csv(sig_res, file.path(TABLE_DIR, "LUAD_rMos_DESeq2_significant.csv"), row.names = FALSE)
write.csv(as.data.frame(pseudobulk_counts) %>% rownames_to_column("gene"),
          file.path(TABLE_DIR, "LUAD_rMos_pseudobulk_counts.csv"), row.names = FALSE)
write.csv(sample_meta %>% rownames_to_column("sample_id"),
          file.path(TABLE_DIR, "LUAD_rMos_pseudobulk_metadata.csv"), row.names = FALSE)

# Variance-stabilized expression is retained for the overlap heatmap.
vsd <- vst(dds, blind = FALSE)
normalized_matrix <- assay(vsd)
saveRDS(normalized_matrix, file.path(OBJECT_DIR, "LUAD_rMos_DESeq2_vst_matrix.rds"))

cat("Significant DESeq2 genes (adjusted p < 0.05): ", nrow(sig_res), "\n", sep = "")
