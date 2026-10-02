# Data

Raw and intermediate data are **not tracked**. Download the files below

and place them under the matching sub-directories before running the

pipeline.

## `data/tcga_kirc/`

| File | Source |

|---|---|

| `KIRC_expr_log2_tpm.csv` | [TCGA GDC](https://portal.gdc.cancer.gov/) project TCGA-KIRC, log2(TPM+1) |

| `KIRC_expr_tpm.csv` | Same source, raw TPM |

## `data/processed/`

| File | Source |

|---|---|

| `subtype_assignment_balanced.csv` | NCC subtype labels (`Patient`, `Subtype`) |

| `deg_limma_voom_with_symbols.csv` | From `Figure3_transcriptomic_CPTAC` |

| `Supplementary_Table_NCC_68_genes.xlsx` | From `Figure1_gene_selection` |

| `cibersort_results_official_full.csv` | CIBERSORT output (spillover-corrected) |

| `ssgsea_immune_scores_R.csv`, `ssgsea_immune_stats_R.csv` | ssGSEA scores and stats from R |

| `immune_genes_stratified_limma_v3.csv` | Stratified limma output |

## `data/xcell/`

| File | Source |

|---|---|

| `xcell_scores_raw.csv` | R xCell v1.1.0 raw scores |

| `xcell_scores_zscore.csv` | Row-wise z-scores |

| `xcell_subtype_diff.csv` | S1 vs S2 differential table |

## `data/scrnaseq/`

GSE159115 ccRCC scRNA-seq:

- `GSE159115_ccRCC_anno.csv.gz` — cell annotations

- `GSM*_SI_*.h5` — one 10x HDF5 file per sample (7 samples)

Download from [GEO GSE159115](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159115).

## `data/gene_sets/`

MSigDB gene sets (symbols):

- `h.all.v2024.1.Hs.symbols.gmt`

- `c2.all.v2024.1.Hs.symbols.gmt`

- `c5.all.v2023.2.Hs.symbols.gmt`

Download from [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/).

## Sibling-repo dependencies

The following two files are produced by other repositories in this

project. Copy them into `data/processed/` before running:

```bash

copy ..\Figure3_transcriptomic_CPTAC\data\processed\deg_limma_voom_with_symbols.csv data\processed\

copy ..\Figure1_gene_selection\results\Supplementary_Table_NCC_68_genes.xlsx data\processed\

```

Alternatively, set the environment variables `LIMMA_CSV` and

`PANEL68_XLSX` to point at those files in place:

```bat

set LIMMA_CSV=D:\github0918\Figure3_transcriptomic_CPTAC\data\processed\deg_limma_voom_with_symbols.csv

set PANEL68_XLSX=D:\github0918\Figure1_gene_selection\results\Supplementary_Table_NCC_68_genes.xlsx

```

The `run_all.py` runner already sets both of these when spawning child

processes.
