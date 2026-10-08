# Copy this file to config/paths.R and edit paths if your local data layout differs.
# The data directory itself is intentionally excluded from version control.

DATA_DIR <- "data"
RESULTS_DIR <- "results"
FIGURE_DIR <- "figures"

# BRCA (SCP1106)
BRCA_COUNTS <- file.path(DATA_DIR, "BRCA", "Wu_EMBO_countr_matrix.csv")
BRCA_METADATA <- file.path(DATA_DIR, "BRCA", "Wu_EMBO_metadata.csv")

# CRC (GSE178341)
CRC_COUNTS_H5 <- file.path(DATA_DIR, "CRC", "GSE178341_crc10x_full_c295v4_submit.h5")
CRC_METADATA <- file.path(DATA_DIR, "CRC", "GSE178341_crc10x_full_c295v4_submit_metatables.csv")

# Melanoma (GSE123139)
MELANOMA_10X_DIR <- file.path(DATA_DIR, "melanoma", "filtered_feature_bc_matrix")
MELANOMA_METADATA <- file.path(DATA_DIR, "melanoma", "Li2018_metadata.csv")

# LUAD (GSE131907)
# Put each sample's 10x matrix directory beneath normal/ or tumor/ using the names
# filtered_feature_bc_matrix.<sample_id>, e.g. filtered_feature_bc_matrix.p018n.
LUAD_NORMAL_DIR <- file.path(DATA_DIR, "LUAD", "normal")
LUAD_TUMOR_DIR <- file.path(DATA_DIR, "LUAD", "tumor")
