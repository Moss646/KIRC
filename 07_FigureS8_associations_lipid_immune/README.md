# 07 — Lipid-metabolism x immune-feature correlation

Computes Spearman, partial Spearman (adjusted for xCell ImmuneScore and
StromaScore), and subtype-stratified Spearman correlations between 30
lipid-metabolism genes and 7 immune features (TIDE + 6 checkpoint genes),
and produces the Figure 6 heatmap + core-association dot plot.

## Pipeline
cd scripts
python run_all.py

text

| Step | Script | Purpose |
|---|---|---|
| 1 | `metab_immune_corr_v4.py` | Compute 210 correlation pairs (30 genes x 7 immune features) across four strata |
| 2 | `fig_Figure6_main_v1.py` | Render Figure 6 (heatmap + core-association dot plot) |

## Inputs required (place them yourself)

Paths are exactly what `config.py` expects.

| File | Place at | Source |
|---|---|---|
| `KIRC_expr_log2_tpm.csv` | `data/tcga_kirc/` | `09_TCGA-KIRC_prep` |
| `subtype_assignment_balanced.csv` | `data/` | `01_Figure1_clustering` |
| `TIDE_official_results.csv` | `data/` | `08_FigureS9_therapeutic_analysis_pubsize` |
| `xcell_scores_raw.csv` | `data/` | `04_Figure4_immune` |

`run_all.py` sets `MEL_DATA_ROOT` to this project's `data/`, so the
expression matrix is resolved as `data/tcga_kirc/KIRC_expr_log2_tpm.csv`.

## Manual copy

```bash
cp ../09_TCGA-KIRC_prep/data/tcga_kirc/KIRC_expr_log2_tpm.csv data/tcga_kirc/
cp ../01_Figure1_clustering/data/subtype_assignment_balanced.csv data/
cp ../08_FigureS9_therapeutic_analysis_pubsize/data/processed/TIDE_official_results.csv data/
cp ../04_Figure4_immune/data/xcell/xcell_scores_raw.csv data/
Outputs
File	Consumed by
results/metabolic_immune_corr_30genes_v1.csv	step 2; downstream meta-analysis
results/Figure6_metab_immune_corr_v1.png	manuscript Figure 6
results/Figure6_metab_immune_corr_v1.svg	manuscript Figure 6
Environment
Python >= 3.10: pandas, numpy, scipy, matplotlib.

Notes
The pipeline asserts that exactly 533 samples remain after merging all
four sources (172 S1, 361 S2). If your upstream files were regenerated
with different filters, this assertion will fail.

fig_Figure6_main_v1.py asserts len(df) == 210 (30 genes x 7 immune
features).

run_all.py sets NOPAUSE=1.