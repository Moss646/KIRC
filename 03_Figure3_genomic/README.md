# 03 — Genomic analysis: mutations, CNA, and FAO9

Processes TCGA-KIRC MAF (mutations) and GISTIC (copy number) files, tests
for subtype associations, and couples CNV with mRNA expression. Produces
Figure 4 and supplementary driver-mutation / 3p-deletion panels.

> **Naming note**: this project is numbered 03 in the repository layout,
> but its scripts are named `fig4_*` because the internal figure is
> "Figure 4" in the manuscript.

## Pipeline
cd scripts
python run_all.py

text

| Step | Script | Purpose |
|---|---|---|
| 1 | `make_fao9_input.py` | Extract TCGA rows from the 9-gene non-NCC FAO score table |
| 2 | `fig4_analysis.py` | Mutation waterfall + CNA frequency (14-gene panel) |
| 3 | `cnv_lipid_remodeling_16genes_check_v1.py` | CNA frequency + CNV-mRNA coupling for the 16 phospholipid-remodeling genes |
| 4 | `driver_mut_fao_v1.py` | Driver mutations vs S1/S2 + FAO9 association + 3p deletion rates |
| 5 | `fig_driver_3p_del_FAO9_v1.py` | Supplementary two-panel figure |
| 6 | `Fig4_plot_v30_pubsize.py` | Figure 4 publication-size (3-panel) |

## Inputs required (place them yourself)

| File | Place at | Source |
|---|---|---|
| `TCGA-KIRC.merged.maf` | `data/raw/` | GDC / FireBrowse download |
| `TCGA-KIRC.gistic.tsv` | `data/raw/` | GDC / Xena download |
| `KIRC_expr_log2_tpm.csv` | `data/raw/` | `09_TCGA-KIRC_prep` |
| `nonncc9_fao_score_per_sample_v1.csv` | `data/raw/` | `05_Figure5_External_validation_OS_KM_score_combined` |
| `subtype_assignment_balanced.csv` | `data/processed/` | `01_Figure1_clustering` |

## Manual copy

```bash
cp ../09_TCGA-KIRC_prep/data/tcga_kirc/KIRC_expr_log2_tpm.csv data/raw/
cp ../01_Figure1_clustering/data/subtype_assignment_balanced.csv data/processed/
cp ../05_Figure5_External_validation_OS_KM_score_combined/results/tables/nonncc9_fao_score_per_sample_v1.csv data/raw/
Outputs
File	Consumed by
data/processed/fao9_immune_input_tcga_v1.csv	driver_mut_fao_v1.py
data/processed/fig4_mutation_stats.csv	Fig4_plot_v30_pubsize.py
data/processed/fig4_oncoprint.csv	Fig4_plot_v30_pubsize.py
data/processed/fig4_cnv_stats.csv	Fig4_plot_v30_pubsize.py
data/processed/fig4_panelA_samples.csv	Fig4_plot_v30_pubsize.py
results/tables/cnv_lipid_remodeling_16genes_check_v2.csv	Figure 4, project 04
results/tables/driver_mut_fao_v1.csv	fig_driver_3p_del_FAO9_v1.py
results/tables/driver_gistic_del_v1.csv	fig_driver_3p_del_FAO9_v1.py
results/figures/Figure_4_genomic_v31_pubsize.png/.svg	manuscript Figure 4
results/figures/Supplementary_Fig_driver_3p_del_FAO9_v1.png/.svg	supplementary
Environment
Python >= 3.10: pandas, numpy, scipy, statsmodels, matplotlib.

Notes
cnv_lipid_remodeling_16genes_check_v1.py calls the Ensembl REST API to
fetch chromosome arms. Without network access the arm column is NA
and the script still runs.

run_all.py sets MEL_DATA_ROOT and NOPAUSE=1.

fig4_params.py must stay in scripts/ next to the other scripts.