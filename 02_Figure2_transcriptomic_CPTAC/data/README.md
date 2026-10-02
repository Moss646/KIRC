# Data

## Files shipped with this repository

| Path | Description |

|---|---|

| `gene_sets/msigdb_lipid_dedup_237.gmt` | 237 deduplicated lipid gene sets (MSigDB-derived), used by Panel C/D/E |

A copy of `msigdb_lipid_dedup_237.gmt` is also expected in

`data/processed/` by `Fig3_main_v18.py`. Place it there too, or symlink it.

## Files NOT shipped with this repository

### `data/raw/`

| File | Size | Source |

|---|---|---|

| `KIRC_expr_log2_tpm.csv` | ~174 MB | [TCGA GDC](https://portal.gdc.cancer.gov/) — project TCGA-KIRC, RNA-seq, log2(TPM+1) |

| `KIRC_counts_for_limma.csv` | large | [TCGA GDC](https://portal.gdc.cancer.gov/) — raw gene counts |

| `rnaseq_tumor.txt` | ~10 MB | [CPTAC KIRC](https://proteomic.datacommons.cancer.gov/pdc/) |

| `proteome_tumor.txt` | ~10 MB | [CPTAC KIRC](https://proteomic.datacommons.cancer.gov/pdc/) |

| `proteome_normal.txt` | ~10 MB | [CPTAC KIRC](https://proteomic.datacommons.cancer.gov/pdc/) |

### `data/processed/`

Produced by the upstream R pipeline in `scripts/upstream/R/`. Rerun it on

the same expression matrix, or request the CSVs from the corresponding

author.

- `subtype_assignment_balanced.csv` — NCC subtype labels

- `subtype_classifier.json` — trained 68-gene NCC classifier

- `balanced_genes.txt` — 68 gene symbols used by the classifier

- `deg_limma_voom.csv` — limma-voom differential expression

- `deg_limma_voom_with_symbols.csv` — DEG with HGNC symbols

- `gsea_full_results.csv` — GSEA results

- `gsva_subtype_statistics.csv` — GSVA pathway statistics

- `gsva_scores.csv` — GSVA per-sample scores

- `msigdb_lipid_dedup_237.gmt` — 237 deduplicated lipid gene sets

## Format requirements

### `KIRC_expr_log2_tpm.csv`

- CSV, comma-separated, UTF-8

- First column: gene symbol; remaining columns: sample barcodes

- Values: log2(TPM+1)

### `proteome_tumor.txt` / `proteome_normal.txt`

- Tab-separated

- First column: gene symbol (index)

- Remaining columns: CPTAC sample IDs

- Values: log2 protein abundance

### `rnaseq_tumor.txt`

- Tab-separated, same orientation as the proteome files

### `deg_limma_voom.csv`

- Must contain columns: `gene`, `logFC`, `AveExpr`, `t`, `p_value`,

&#x20; `FDR`, `category`

- `category` is one of `UP_in_S1`, `DOWN_in_S1`, `NS`

### `deg_limma_voom_with_symbols.csv`

- Same columns as `deg_limma_voom.csv` plus `gene_symbol`

### `subtype_assignment_balanced.csv`

- Columns: `Patient`, `Subtype` (`Subtype_1` or `Subtype_2`)

### `subtype_classifier.json`

- Keys: `genes`, `centroid_s1`, `centroid_s2`, `scaler_mean`, `scaler_std`
