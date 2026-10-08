# Matched-patient LUAD monocytic subset abundance analysis.
# Relative abundance is calculated within each patient x tissue specimen.

source(file.path("scripts", "00_setup.R"))
mono <- readRDS(file.path(OBJECT_DIR, "LUAD_high_resolution_monocytes.rds"))

subset_levels <- levels(factor(mono$mono_subsets))
observed_samples <- mono@meta.data %>%
  distinct(patient_id, tissue_type)
skeleton <- observed_samples %>%
  tidyr::crossing(mono_subsets = subset_levels)

abundance <- mono@meta.data %>%
  count(patient_id, tissue_type, mono_subsets, name = "n_cells") %>%
  right_join(skeleton, by = c("patient_id", "tissue_type", "mono_subsets")) %>%
  mutate(n_cells = replace_na(n_cells, 0L)) %>%
  group_by(patient_id, tissue_type) %>%
  mutate(total_monocytes = sum(n_cells),
         relative_abundance = if_else(total_monocytes > 0, n_cells / total_monocytes, NA_real_)) %>%
  ungroup()

# Keep patients represented in both Normal and Tumor samples.
paired_ids <- abundance %>%
  distinct(patient_id, tissue_type) %>%
  count(patient_id) %>%
  filter(n == 2) %>%
  pull(patient_id)
paired <- abundance %>% filter(patient_id %in% paired_ids)

paired_wide <- paired %>%
  select(patient_id, tissue_type, mono_subsets, relative_abundance) %>%
  pivot_wider(names_from = tissue_type, values_from = relative_abundance)

stats <- paired_wide %>%
  group_by(mono_subsets) %>%
  summarise(
    n_pairs = sum(complete.cases(Normal, Tumor)),
    normal_mean = mean(Normal, na.rm = TRUE),
    tumour_mean = mean(Tumor, na.rm = TRUE),
    p_value = if (sum(complete.cases(Normal, Tumor)) > 0) {
      wilcox.test(Tumor, Normal, paired = TRUE, exact = FALSE)$p.value
    } else NA_real_,
    .groups = "drop"
  )

write.csv(abundance, file.path(TABLE_DIR, "LUAD_monocytic_subset_relative_abundance.csv"), row.names = FALSE)
write.csv(stats, file.path(TABLE_DIR, "LUAD_monocytic_subset_paired_Wilcoxon.csv"), row.names = FALSE)

p <- ggplot(paired, aes(x = mono_subsets, y = relative_abundance, group = interaction(patient_id, tissue_type))) +
  geom_boxplot(aes(group = interaction(mono_subsets, tissue_type)), outlier.shape = NA) +
  geom_jitter(aes(shape = tissue_type), width = 0.12, height = 0, alpha = 0.8) +
  labs(x = NULL, y = "Relative abundance") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p, "LUAD_monocytic_subset_abundance.pdf", width = 10, height = 6)

print(stats)
