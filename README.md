# TCGA-KIRC Lipid-Metabolism Subtype Multi-Omics Analysis

This repository contains nine projects that together build and validate a
lipid-metabolism subtype classifier for TCGA-KIRC clear cell renal cell
carcinoma, with transcriptomic, proteomic, genomic, immune, therapeutic,
and external cohort validation.

## Repository layout

| Folder | Content |
|---|---|
| `01_Figure1_clustering/` | Consensus clustering, 68-gene nearest-centroid classifier (NCC), NMF validation, survival-blind sensitivity analysis |
| `02_Figure2_transcriptomic_CPTAC/` | limma-voom differential expression, GSEA/GSVA, CPTAC proteomic validation |
| `03_Figure3_genomic/` | Driver mutations, GISTIC copy-number alterations, lipid-gene CNA checks |
| `04_Figure4_immune/` | Immune deconvolution (CIBERSORT, ssGSEA, xCell), scRNA-seq (GSE159115), HIF/hypoxia analyses |
| `05_Figure5_External_validation_OS_KM_score_combined/` | External validation in E-MTAB-1980, ICGC RECA-EU, CPTAC: Cox models, KM curves, cross-cohort concordance, tables |
| `06_FigureS1_gene_selection/` | Lipid gene-set collection (MSigDB + GO dedup), Cox screening, 68-gene panel |
| `07_FigureS8_associations_lipid_immune/` | Correlation between lipid-metabolism genes and immune features (xCell, TIDE) |
| `08_FigureS9_therapeutic_analysis_pubsize/` | oncoPredict/GDSC2 drug-sensitivity prediction and TIDE immune-response scores |
| `09_TCGA-KIRC_prep/` | Raw TCGA-KIRC download and preprocessing (expression, clinical, counts) |

## Dependency order
09_TCGA-KIRC_prep (tcga_kirc_prep.py: expression / clinical)
|
v
06_FigureS1_gene_selection (lipid gene universe -> 68-gene panel)
|
v
01_Figure1_clustering (S1/S2 subtypes + frozen subtype_classifier.json)
|
v
09_TCGA-KIRC_prep (download_gdc_counts.py: STAR counts for limma)
|
+--> 02_Figure2_transcriptomic_CPTAC --> (deg table consumed by 04)
+--> 03_Figure3_genomic
+--> 04_Figure4_immune
+--> 05_Figure5_External_validation_OS_KM_score_combined
+--> 07_FigureS8_associations_lipid_immune
+--> 08_FigureS9_therapeutic_analysis_pubsize

text

## Environment

- **Python** >= 3.10 (tested with 3.13). See `requirements.txt`.
- **R** >= 4.5.0 (tested with 4.5.0). Packages used: `limma`, `edgeR`,
  `GSVA`, `fgsea`, `org.Hs.eg.db` (Bioconductor); `xCell` (GitHub);
  `CIBERSORT` (registration required); `oncoPredict` (Bioconductor/GitHub).
- Most Python scripts resolve their paths relative to their own location
  and need no environment variables. Projects 03, 04, 05, and 07 set or
  honor `MEL_DATA_ROOT` to point at their own `data/` directory; projects
  03 and 04 also set `NOPAUSE=1` to skip the interactive "press any key"
  prompt when run non-interactively.

## Data policy

Large raw/intermediate data files are **not tracked** in this repository
(GitHub 100 MB per-file limit). Each project README lists every required
input file, the exact repo-relative location it must be placed at, and
the public source to obtain it. Place the files first, then run the
pipelines.

## Running a project

Every project contains a `scripts/run_all.py` (or `scripts/run_upstream.py`)
that executes the full pipeline in the correct order:
cd 01_Figure1_clustering/scripts
python run_all.py

text

Outputs are written to `results/` (figures, tables) and, where applicable,
to `data/processed/` (intermediate files consumed by later projects).

## Note on duplicated filenames

`driver_mut_fao_v1.py` exists in both `03_Figure3_genomic/` and
`05_Figure5_External_validation_OS_KM_score_combined/`. The two versions
are functionally equivalent; project 05 is the canonical one, and project
03 uses its own copy only for producing the supplementary figure. Do not
confuse them.

## License

MIT — see `LICENSE`.
