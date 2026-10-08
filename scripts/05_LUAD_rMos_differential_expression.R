# Single-cell differential expression between LUAD tumour-associated and normal-tissue rMos.

source(file.path("scripts", "00_setup.R"))
mono <- readRDS(file.path(OBJECT_DIR, "LUAD_high_resolution_monocytes.rds"))

mono$rMos_condition <- paste(mono$mono_subsets, mono$tissue_type, sep = "_")
Idents(mono) <- "rMos_condition"

rmos_de <- FindMarkers(
  mono,
  ident.1 = "rMos_Tumor",
  ident.2 = "rMos_Normal",
  test.use = "wilcox",
  verbose = FALSE
)
rmos_de <- rmos_de %>% rownames_to_column("gene")
rmos_de$significant <- rmos_de$p_val_adj < 0.05
rmos_de$volcano_highlight <- rmos_de$p_val_adj < 0.05 & abs(rmos_de$avg_log2FC) > 2

write.csv(rmos_de, file.path(TABLE_DIR, "LUAD_rMos_FindMarkers_all.csv"), row.names = FALSE)
write.csv(filter(rmos_de, significant),
          file.path(TABLE_DIR, "LUAD_rMos_FindMarkers_significant.csv"), row.names = FALSE)

# Volcano plot uses adjusted p-values on the y-axis.
plot_df <- rmos_de %>%
  mutate(log10_adj_p = -log10(pmax(p_val_adj, .Machine$double.xmin)),
         class = case_when(
           volcano_highlight & avg_log2FC > 0 ~ "Up",
           volcano_highlight & avg_log2FC < 0 ~ "Down",
           TRUE ~ "Not highlighted"
         ))

p_volcano <- ggplot(plot_df, aes(avg_log2FC, log10_adj_p)) +
  geom_point(aes(shape = class), alpha = 0.65, size = 1.4) +
  geom_vline(xintercept = c(-2, 2), linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  labs(x = "Average log2 fold change", y = expression(-log[10]("adjusted p-value")),
       title = "LUAD rMos: Tumor vs Normal") +
  theme_classic()
save_plot(p_volcano, "LUAD_rMos_volcano.pdf", width = 7, height = 6)

cat("Significant Seurat DEGs (adjusted p < 0.05): ", sum(rmos_de$significant, na.rm = TRUE), "\n", sep = "")
