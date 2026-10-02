# 08 — Therapeutic analysis: drug sensitivity and TIDE

Predicts per-sample drug sensitivity (IC50) for nine targeted agents via
oncoPredict::calcPhenotype (GDSC2 training), runs TIDE via tidepy to score
predicted immune response, and produces Figure 7 plus a faceted per-drug
panel.

> **Naming note**: this project is numbered 08/FigureS9 in the repository
> layout, but its scripts are named `Fig7_*` because the internal figure
> is "Figure 7" in the manuscript.

## Pipeline
cd scripts
python run_all.py

text

| Step | Script | Purpose |
|---|---|---|
| 1 | `prepare_oncopredict_inputs.py` | Copy `KIRC_expr_log2_tpm.csv` and `subtype_assignment_balanced.csv` into `data/raw/` as `kirc_expr_matched.csv` and `kirc_subtypes.csv` |
| 2 | `run_oncoPredict_calcPheno.R` | Predict IC50 for 9 drugs via oncoPredict |
| 3 | `run_tidepy.py` | Run tidepy on the KIRC expression matrix |
| 4 | `Fig7_main_pubsize_v2.py` | Render Figure 7 (six panels) |
| 5 | `Fig7_PanelA_violin_faceted.py` | Faceted per-drug violin/box panel |

## Inputs required (place them yourself)

| File | Place at | Source |
|---|---|---|
| `KIRC_expr_log2_tpm.csv` | `data/` or `data/raw/` | `09_TCGA-KIRC_prep` |
| `subtype_assignment_balanced.csv` | `data/` or `data/raw/` | `01_Figure1_clustering` |
| `GDSC2_Expr.rds` | `data/raw/` | GDSC2 / oncoPredict documentation |
| `GDSC2_Res.rds` | `data/raw/` | GDSC2 / oncoPredict documentation |

## Manual copy

```bash
cp ../09_TCGA-KIRC_prep/data/tcga_kirc/KIRC_expr_log2_tpm.csv data/
cp ../01_Figure1_clustering/data/subtype_assignment_balanced.csv data/
prepare_oncopredict_inputs.py scans both data/ and data/raw/ and
writes canonical copies to data/raw/kirc_expr_matched.csv and
data/raw/kirc_subtypes.csv.

Outputs
File	Consumed by
data/processed/oncoPredict_calcPheno_results.csv	step 4
data/processed/oncoPredict_calcPheno_ic50.csv	steps 4, 5
data/processed/TIDE_official_results.csv	step 4; 07_FigureS8_associations_lipid_immune
results/figures/Figure_7_therapeutic_analysis_pubsize_v2.png/.svg	manuscript Figure 7
results/figures/Fig7_PanelA_violin_faceted.png/.svg	per-drug faceted panel
Environment
Python >= 3.10: pandas, numpy, scipy, statsmodels, matplotlib,
Pillow, plus tidepy (GitHub, not PyPI).
R >= 4.5.0: oncoPredict, glmnet, ridge.
Rscript must be on PATH.

Notes
run_all.py sets MEL_DATA_ROOT to the project root and runs all five
steps via subprocess.

oncoPredict writes a temporary calcPhenotype_Output/ folder inside
data/processed/; the R script removes it at the end.

Drug order in Panels A/B/C is determined by ascending residual Cohen's
d.