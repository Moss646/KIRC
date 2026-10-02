# 01 — Consensus clustering and the 68-gene NCC classifier

Runs consensus clustering of the TCGA-KIRC cohort on the 68-gene lipid panel
(two subtypes, S1/S2), trains the frozen 68-gene nearest-centroid classifier
(NCC), validates it with NMF, performs the survival-blind sensitivity
analysis, and produces the clustering figure, baseline table, and classifier
table.

## Pipeline
cd scripts
python run_upstream.py # steps 1-13 below, in order
python run_upstream.py --list # inspect / step-select

text

| Step | Script | Purpose |
|---|---|---|
| 1 | `compute_clustering.py` | Consensus clustering (writes `cache/consensus_cache.pkl` and `data/subtype_assignment_balanced.csv`) |
| 2 | `build_classifier.py` | Trains `data/subtype_classifier.json`; writes an ICGC validation R skeleton to `data/icgc_validation.R` |
| 3-4 | `plot_figure2.py`, `plot_figure2_pubsize.py` | Main clustering figure (standard + publication size) |
| 5-6 | `supp_cdf_plot.py`, `supp_cdf_plot_pubsize.py` | CDF / consensus-index supplementary plots |
| 7-9 | `nmf_validation.py`, `nmf_validation_metrics_v2.py`, `nmf_validation_pubsize.py` | NMF (vs k-means) validation |
| 10-11 | `table1_baseline.py`, `table_s2_ncc_classifier.py` | Baseline characteristics and classifier tables |
| 12 | `survival_blind_analysis_v1.py` | Survival-blind sensitivity analysis (~9 min; k = 2-6 consensus clustering) |
| 13 | `survival_blind_supplementary_figure_v4.py` | 4-panel supplementary figure |

## Inputs required (place them yourself)

The `data/` directory ships **empty** — every input must be placed before
running.

| File | Place at | Source | Used by |
|---|---|---|---|
| `KIRC_expr_log2_tpm.csv` | `data/` | `09_TCGA-KIRC_prep` output (`data/tcga_kirc/`) | all steps |
| `balanced_genes.txt` | `data/` | `06_FigureS1_gene_selection` output (`data/consensus_cluster/`) | all steps |
| `TCGA-KIRC_clinical_info.csv` | `data/` | GDC clinical export | steps 1, 10 |
| `TCGA-KIRC_clinical_survival.csv` | `data/` | GDC clinical export | steps 1, 2 |
| `KIRC_clinical_cleaned.csv` | `data/` | `09_TCGA-KIRC_prep` output | step 12 |
| `msigdb_lipid_gene_pool_dedup.csv` | `data/` | `06_FigureS1_gene_selection` output | step 12 |

`subtype_assignment_balanced.csv` is **generated** here by step 1 on first
run.

## Manual copy

Before running, copy upstream outputs into `data/`.

```bash
cp ../09_TCGA-KIRC_prep/data/tcga_kirc/KIRC_expr_log2_tpm.csv       data/
cp ../09_TCGA-KIRC_prep/data/tcga_kirc/KIRC_clinical_cleaned.csv    data/
cp ../06_FigureS1_gene_selection/data/consensus_cluster/balanced_genes.txt  data/
cp ../06_FigureS1_gene_selection/data/msigdb_lipid_gene_pool_dedup.csv      data/
# raw GDC clinical exports (user-supplied):
#   data/TCGA-KIRC_clinical_info.csv
#   data/TCGA-KIRC_clinical_survival.csv
Outputs
File	Consumed by
data/subtype_assignment_balanced.csv	02, 03, 04, 05, 07, 08
data/subtype_classifier.json	02, 04, 05
data/icgc_validation.R	optional ICGC external validation
results/ figures and tables	manuscript Figure 2, CDF supplement, NMF validation, Table 1, Table S2
Optional: ICGC external validation
build_classifier.py writes data/icgc_validation.R. To use it, download
donor.RECA-EU.tsv.gz, specimen.RECA-EU.tsv.gz, and
exp_seq.RECA-EU.tsv.gz from
ICGC, put them in
the same directory as the R script together with subtype_classifier.json,
and run it with R.

Environment
Python >= 3.10: pandas, numpy, scipy, scikit-learn, lifelines,
matplotlib, statsmodels. Fixed random seed (20261002) for the
survival-blind analysis; consensus clustering caches results in cache/.

Note
subtype_assignment_balanced.csv and subtype_classifier.json are the two
key artifacts that all downstream projects consume. If you rerun this
project and the cached files change, downstream projects must be rerun with
the newly produced copies.