# Data

## Bundled in this repository

| File | Description |
|---|---|
| `subtype_assignment_balanced.csv` | Columns: `Patient`, `Subtype` (`Subtype_1` or `Subtype_2`) |
| `TIDE_official_results.csv` | TIDE scores per sample (first column = sample ID) |
| `xcell_scores_raw.csv` | xCell cell-type enrichment (first column = `celltype`) |

## Not bundled (required externally)

| File | Source |
|---|---|
| `KIRC_expr_log2_tpm.csv` (~174 MB) | [TCGA GDC](https://portal.gdc.cancer.gov/) project TCGA-KIRC, log2(TPM+1) |

Place it at `<MEL_DATA_ROOT>/tcga_kirc/KIRC_expr_log2_tpm.csv` and set
`MEL_DATA_ROOT` to `<MEL_DATA_ROOT>` before running the scripts. See the
top-level `README.md` for the exact commands.
