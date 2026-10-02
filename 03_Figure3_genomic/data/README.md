# Data

Raw and intermediate data are **not tracked**. Place the files below
under the matching sub-directories.

## `data/raw/`

| File | Source |
|---|---|
| `TCGA-KIRC.merged.maf` | [TCGA GDC](https://portal.gdc.cancer.gov/) — project TCGA-KIRC, somatic mutations |
| `TCGA-KIRC.gistic.tsv` | [TCGA GDC](https://portal.gdc.cancer.gov/) / [UCSC Xena](https://xenabrowser.net/) — GISTIC2.0 gene-level CNA |
| `KIRC_expr_log2_tpm.csv` | [TCGA GDC](https://portal.gdc.cancer.gov/) — gene × sample log2(TPM+1) |

### Formats

- **`TCGA-KIRC.merged.maf`** — tab-separated MAF. Columns used:
  `Tumor_Sample_Barcode`, `Hugo_Symbol`, `Variant_Classification`.
  Comment lines start with `#`.
- **`TCGA-KIRC.gistic.tsv`** — gene × sample matrix. First column is the
  gene symbol; remaining columns are TCGA sample barcodes.
- **`KIRC_expr_log2_tpm.csv`** — gene × sample matrix. First column is
  the gene symbol; remaining columns are TCGA sample barcodes.

## `data/processed/`

| File | Source |
|---|---|
| `subtype_assignment_balanced.csv` | Upstream NCC classifier (see `Figure2_clustering`) |
| `fao9_immune_input_tcga_v1.csv` | Upstream FAO9 scoring (nine-gene score per sample) |

`subtype_assignment_balanced.csv` must contain at least:
- `Patient` — 12-character TCGA patient ID
- `Subtype` — one of `Subtype_1`, `Subtype_2`

`fao9_immune_input_tcga_v1.csv` must contain at least:
- `sample` — 12-character TCGA patient ID
- `FAO9` — the nine-gene FAO score

## Overriding the expression-matrix path

If `KIRC_expr_log2_tpm.csv` is stored outside the repository, set the
`KIRC_EXPR_CSV` environment variable to its absolute path:

```bash
# Windows CMD
set KIRC_EXPR_CSV=D:\path\to\KIRC_expr_log2_tpm.csv
# Linux / macOS
export KIRC_EXPR_CSV=/path/to/KIRC_expr_log2_tpm.csv
```
