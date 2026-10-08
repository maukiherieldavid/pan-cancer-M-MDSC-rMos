# Pan-cancer M-MDSC signature and rMos analysis

This repository contains the cleaned R analysis workflow associated with the manuscript:

**A conserved pan-cancer transcriptional signature supports identification of human M-MDSCs and reveals REL/NFKB1⁺ rMos across human cancers**

The repository is a reproducible methodological representation of the analyses reported in the manuscript. Public source datasets are not redistributed.

## Study overview

The workflow uses breast cancer single-cell RNA-sequencing data for discovery of a candidate five-gene M-MDSC signature (`VCAN`, `APOBEC3A`, `CD300E`, `EREG`, `FCN1`) and independently evaluates the signature in colorectal cancer, melanoma, and lung adenocarcinoma. Higher-resolution monocytic analyses are then used to identify REL/CCR2-associated rMo-like populations. LUAD rMos are further analysed by single-cell differential expression, sample-level pseudobulk DESeq2, overlap analysis, pathway/GO enrichment, and matched-patient monocytic subset abundance analysis.

## Public datasets

| Cancer type | Accession | Study role |
|---|---|---|
| Breast cancer (BRCA) | SCP1106 | Signature discovery |
| Colorectal cancer (CRC) | GSE178341 | Independent transcriptomic evaluation |
| Melanoma | GSE123139 | Independent transcriptomic evaluation |
| Lung adenocarcinoma (LUAD) | GSE131907 | Independent transcriptomic evaluation and rMos analysis |

The melanoma data are also distributed by the source study at `https://github.com/tanaylab/li_et_al_cell_2018_melanoma_scrna`.

## Software

The manuscript analysis used **R 4.3.0** and **Seurat 5.1.0**. Additional packages used by the cleaned workflow include `dplyr`, `tidyr`, `tibble`, `stringr`, `readr`, `ggplot2`, `Matrix`, `patchwork`, `DESeq2`, `apeglm`, `clusterProfiler`, `org.Hs.eg.db`, `ComplexHeatmap`, and `circlize`.

Record the local environment after completing the analysis with:

```r
source("scripts/99_session_info.R")
```

## Repository structure

```text
.
├── README.md
├── config/
│   └── paths.example.R
├── metadata/
│   ├── datasets.tsv
│   ├── expected_results.tsv
│   ├── external_analyses.md
│   └── README.md
├── scripts/
│   ├── 00_setup.R
│   ├── 01_BRCA_signature_discovery.R
│   ├── 02_CRC_signature_evaluation.R
│   ├── 03_melanoma_signature_evaluation.R
│   ├── 04_LUAD_signature_and_rMos.R
│   ├── 05_LUAD_rMos_differential_expression.R
│   ├── 06_LUAD_pseudobulk_DESeq2.R
│   ├── 07_LUAD_overlap_enrichment.R
│   ├── 08_LUAD_monocytic_subset_abundance.R
│   ├── 09_computational_figure_panels.R
│   └── 99_session_info.R
├── results/
│   └── README.md
└── figures/
    └── README.md
```

## Setup

1. Download the public datasets from the repositories above.
2. Copy `config/paths.example.R` to `config/paths.R`.
3. Edit `config/paths.R` if your local data layout differs.
4. Run the scripts from the repository root in numerical order.

Example:

```r
source("scripts/01_BRCA_signature_discovery.R")
source("scripts/02_CRC_signature_evaluation.R")
source("scripts/03_melanoma_signature_evaluation.R")
source("scripts/04_LUAD_signature_and_rMos.R")
source("scripts/05_LUAD_rMos_differential_expression.R")
source("scripts/06_LUAD_pseudobulk_DESeq2.R")
source("scripts/07_LUAD_overlap_enrichment.R")
source("scripts/08_LUAD_monocytic_subset_abundance.R")
source("scripts/09_computational_figure_panels.R")
source("scripts/99_session_info.R")
```

## Analysis logic

### BRCA signature discovery

`01_BRCA_signature_discovery.R` uses RPCA integration and resolves the immune-cell compartment. Candidate transcriptional markers are assessed by three complementary analyses:

- `FindAllMarkers` across annotated immune-cell populations (`only.pos = TRUE`, `min.pct = 0.25`, `logfc.threshold = 0.25`).
- `FindMarkers` comparing the two monocytic clusters with macrophages.
- `FindMarkers` comparing the two monocytic clusters with all other annotated immune-cell populations.

For the two targeted comparisons, RPS/RPL ribosomal genes are excluded for visualization and the 25 positive genes with the highest `avg_log2FC` are exported. Candidate genes are then evaluated across the full annotated immune compartment. The final five genes are the candidates showing preferential expression in both monocytic clusters with comparatively limited expression in other immune populations; the five-gene set is therefore **not treated as the automatic intersection of two DEG lists**.

### Cross-cancer evaluation

The CRC, melanoma, and LUAD scripts evaluate `VCAN`, `APOBEC3A`, `CD300E`, `EREG`, and `FCN1` within monocytic populations together with conventional/contextual features including `CD14`, `FCN1`, HLA-DR genes, and `TGFB1`. Higher-resolution analysis is performed on the candidate monocytic compartment, followed by evaluation of REL/CCR2-associated rMo-like populations and NFKB1 expression.

### LUAD rMos differential expression

`05_LUAD_rMos_differential_expression.R` compares tumour-associated and normal-tissue rMos using Seurat `FindMarkers` with the Wilcoxon rank-sum test. Genes with Benjamini-Hochberg-adjusted `p < 0.05` define the significant Seurat DEG set. The volcano plot highlights genes with `|avg_log2FC| > 2` and adjusted `p < 0.05`; its y-axis is `-log10(adjusted p-value)`.

### Pseudobulk DESeq2 and overlapping DEGs

`06_LUAD_pseudobulk_DESeq2.R` aggregates raw rMos counts by sample and performs DESeq2 with tissue type as the design factor. DESeq2 significance is defined as adjusted `p < 0.05`.

`07_LUAD_overlap_enrichment.R` computes the overlap between the significant Seurat and DESeq2 gene sets programmatically. Downstream directionality is assigned using the Seurat `avg_log2FC`. Overlapping genes are filtered at `|Seurat avg_log2FC| > 0.5` before KEGG and GO enrichment analyses.

The reported analysis yields:

- 1,940 significant Seurat DEGs
- 1,336 significant DESeq2 DEGs
- 1,052 genes shared between the two analyses
- 2,224 genes in the union (47.3% overlap)
- 856 filtered overlapping genes for enrichment: 510 upregulated and 346 downregulated

### Matched LUAD monocytic subset abundance

`08_LUAD_monocytic_subset_abundance.R` calculates the abundance of each monocytic subset separately within each patient × tissue specimen, explicitly includes zero-count subsets, restricts the comparison to patients represented in both tumour and normal tissue, and applies a two-sided paired Wilcoxon signed-rank test.

## Notes on generated outputs

Large Seurat objects and source datasets are excluded from version control. Tables and computational figure panels are written to `results/` and `figures/`. The multiplex immunohistochemistry figures are based on microscopy images and are not generated by the R workflow. GEPIA2/OncoDB analyses are documented in `metadata/external_analyses.md`.

## Cluster labels and biological interpretation

Cluster numbers are analysis-specific identifiers. The scripts preserve the cluster identities used in the reported analysis, but biological interpretation is based on marker expression rather than on the numerical label itself. If a source dataset or software version is changed, re-check canonical lineage markers and the candidate M-MDSC/rMo-associated features before transferring cluster annotations.
