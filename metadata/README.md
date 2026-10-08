# Input data

Public single-cell datasets are not redistributed in this repository. Download them from their original repositories and arrange them under `data/` as described in `config/paths.example.R`.

| Cancer type | Accession | Role in study |
|---|---|---|
| Breast cancer | SCP1106 | Five-gene signature discovery |
| Colorectal cancer | GSE178341 | Independent transcriptomic evaluation |
| Melanoma | GSE123139 | Independent transcriptomic evaluation |
| Lung adenocarcinoma | GSE131907 | Independent transcriptomic evaluation and rMos analysis |

The analysis uses the source-study metadata only where required for sample, patient, tissue, or annotation information. For LUAD, `sample_id`, `patient_id`, and `tissue_type` are derived programmatically from the sample directory names.
