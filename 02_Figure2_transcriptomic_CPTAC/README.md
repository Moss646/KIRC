# 02 — Transcriptomic DEG, GSEA/GSVA, and CPTAC validation

Runs limma-voom differential expression (S1 vs S2), GSEA via fgsea, GSVA
pathway scoring, and CPTAC protein-level validation. Produces Figure 3.

## Pipeline
cd scripts
python run_all.py

text

| Step | Script | Purpose |
|---|---|---|
| 1 | `02_differential_exp.R` | edgeR + limma-voom S1 vs S2 DEG |
| 2 | `01_convert_ensembl_to_symbol.R` | Map Ensembl IDs to gene symbols via org.Hs.eg.db |
| 3 | `make_msigdb_gmt.py` | Build the 236-set MSigDB lipid GMT |
| 4 | `03_fgsea.R` | fgsea on the ranked logFC list |
| 5 | `gsva_R.R` | GSVA scores (Gaussian KCDF) |
| 6 | `gsva_limma.R` | limma eBayes on the GSVA scores |
| 7 | `Fig3_main_v18.py` | Figure 3 (3x3 panel layout) |
| 8 | `cptac_lipid_remodeling_check_v2.py` | CPTAC protein check for 16 phospholipid-remodeling genes |

## Inputs required (place them yourself)

| File | Place at | Source |
|---|---|---|
| `KIRC_counts_for_limma.csv` | `data/raw/` | `09_TCGA-KIRC_prep` (download_gdc_counts.py) |
| `sample_info_for_limma.csv` | `data/raw/` | `09_TCGA-KIRC_prep` (download_gdc_counts.py) |
| `KIRC_expr_log2_tpm.csv` | `data/raw/` | `09_TCGA-KIRC_prep` (tcga_kirc_prep.py) |
| `rnaseq_tumor.txt`, `proteome_tumor.txt`, `proteome_normal.txt` | `data/raw/` | CPTAC downloads |
| `msigdb_lipid_genesets_dedup.csv` | `data/gene_sets/` | `06_FigureS1_gene_selection` |
| `subtype_assignment_balanced.csv` | `data/processed/` | `01_Figure1_clustering` |
| `subtype_classifier.json` | `data/processed/` | `01_Figure1_clustering` |

## Manual copy

```bash
cp ../09_TCGA-KIRC_prep/data/tcga_kirc/KIRC_counts_for_limma.csv data/raw/
cp ../09_TCGA-KIRC_prep/data/tcga_kirc/sample_info_for_limma.csv data/raw/
cp ../09_TCGA-KIRC_prep/data/tcga_kirc/KIRC_expr_log2_tpm.csv    data/raw/
cp ../01_Figure1_clustering/data/subtype_assignment_balanced.csv data/processed/
cp ../01_Figure1_clustering/data/subtype_classifier.json         data/processed/
cp ../06_FigureS1_gene_selection/data/msigdb_lipid_genesets_dedup.csv data/gene_sets/
Outputs
File	Consumed by
data/processed/deg_limma_voom.csv	Figure 3
data/processed/deg_limma_voom_with_symbols.csv	04, 05, Figure 3
data/processed/gsea_full_results.csv	Figure 3
data/processed/gsva_scores.csv	04, Figure 3
data/processed/gsva_subtype_statistics.csv	Figure 3
data/gene_sets/msigdb_lipid_dedup_236.gmt	03_fgsea.R, gsva_R.R
results/figures/Figure_3_transcriptomic_CPTAC_v18.png/.svg	manuscript Figure 3
results/tables/cptac_lipid_remodeling_16genes_check_v2.csv	Figure 3, project 03
Environment
Python >= 3.10: pandas, numpy, scipy, statsmodels, scikit-learn,
matplotlib.
R >= 4.5.0: edgeR, limma, fgsea, GSVA, org.Hs.eg.db, AnnotationDbi.
Rscript must be on PATH.

Notes
run_all.py sets MEL_DATA_ROOT to the project root, so all scripts see
data/raw/, data/processed/, data/gene_sets/ under this project.

All R scripts use .. from scripts/ to reach the project root, so they
must be run from inside scripts/ (or via run_all.py).