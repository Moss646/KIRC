# 04 — Immune microenvironment and scRNA-seq validation

Runs xCell, CIBERSORT, ssGSEA (Hallmark immune), three-method HIF
scoring (ssGSEA, GSVA, camera), scRNA-seq pseudobulk subtyping, purity
stratification, and produces Figure 4/5 and the scRNA-seq supplementary
figure.

## Pipeline
cd scripts
python run_all.py

text

The 26 steps are listed in `scripts/run_all.py`. Broadly:

1. xCell spillover extraction + xCell analysis
2. CIBERSORT (LM22, perm=1000)
3. ssGSEA immune + HIF panels (ssGSEA, GSVA, camera)
4. Purity strata + stratified immune MWU
5. scRNA-seq pseudobulk + per-cell-type checkpoint/MHC expression
6. Immune-feature composition independence + lipid remodeling + PLA2 panels
7. Publication-size figures

## Inputs required (place them yourself)

| File | Place at | Source |
|---|---|---|
| `KIRC_expr_log2_tpm.csv` | `data/tcga_kirc/` | `09_TCGA-KIRC_prep` |
| `KIRC_expr_tpm.csv` | `data/tcga_kirc/` | `09_TCGA-KIRC_prep` |
| `subtype_assignment_balanced.csv` | `data/processed/` | `01_Figure1_clustering` |
| `subtype_classifier.json` | `data/processed/` | `01_Figure1_clustering` |
| `deg_limma_voom_with_symbols.csv` | `data/processed/` | `02_Figure2_transcriptomic_CPTAC` |
| `gsva_scores.csv` | `data/processed/` | `02_Figure2_transcriptomic_CPTAC` |
| `selected_genes_strict.csv` | `data/processed/` | `06_FigureS1_gene_selection` |
| MSigDB GMTs | `data/gene_sets/` | Download from MSigDB (h.all, c2.all, c5.all) |
| GSE159115 scRNA-seq | `data/scrnaseq/` | GEO GSE159115 download |

## Manual copy

```bash
cp ../09_TCGA-KIRC_prep/data/tcga_kirc/KIRC_expr_log2_tpm.csv data/tcga_kirc/
cp ../09_TCGA-KIRC_prep/data/tcga_kirc/KIRC_expr_tpm.csv      data/tcga_kirc/
cp ../01_Figure1_clustering/data/subtype_assignment_balanced.csv data/processed/
cp ../01_Figure1_clustering/data/subtype_classifier.json         data/processed/
cp ../02_Figure2_transcriptomic_CPTAC/data/processed/deg_limma_voom_with_symbols.csv data/processed/
cp ../02_Figure2_transcriptomic_CPTAC/data/processed/gsva_scores.csv                 data/processed/
cp ../06_FigureS1_gene_selection/data/gene_selection/selected_genes_strict.csv       data/processed/
Outputs (selected)
File	Consumed by
data/xcell/xcell_scores_raw.csv	07
data/processed/cibersort_results_official_full.csv	Figure 4
data/processed/ssgsea_immune_scores_R.csv	Figure 4/5, 07 (HIF)
results/figures/fig4_immune_merged_v17.png/.svg	manuscript Figure 4
results/figures/Fig_supp_scrnaseq_v5_pubsize.png/.svg	scRNA-seq supplementary
data/processed/Supplementary_Table_NCC_68_genes.xlsx	04 internal (panel68 binning)
Full list of outputs appears under results/tables/ and results/figures/.

Environment
Python >= 3.10: pandas, numpy, scipy, statsmodels, matplotlib,
h5py, openpyxl.
R >= 4.5.0: xCell, CIBERSORT (registration required), GSVA, limma,
edgeR.

Notes
run_all.py sets MEL_DATA_ROOT to the project root. All downstream
scripts read data/raw/, data/processed/, data/gene_sets/,
data/tcga_kirc/, data/xcell/, data/scrnaseq/ under this project.

test_xcell.R is a diagnostic script; it is not part of the main
pipeline.

GSE159115 download is large; grab it separately from GEO and place the
.h5 files and _anno.csv.gz under data/scrnaseq/.