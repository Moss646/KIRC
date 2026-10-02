# 09 — TCGA-KIRC raw data preprocessing

Converts raw TCGA-KIRC downloads (gene annotation, mRNA TPM, clinical) into
the standardized expression matrix and cleaned clinical table consumed by
every downstream project. Also downloads GDC STAR counts for the limma-voom
pipeline.

This project has **two entry scripts**, run at different points in the
overall pipeline:

- `tcga_kirc_prep.py` — runs first (raw -> cleaned)
- `download_gdc_counts.py` — runs after `01_Figure1_clustering`

## Pipeline

### Step A: raw -> cleaned
cd scripts
python tcga_kirc_prep.py

text

Four cleaning steps:

1. keep protein-coding genes only
2. keep primary-tumour samples only (sample code `01`)
3. keep the earliest vial per patient (01A before 01B)
4. low-expression filter: mean TPM >= 1.0 AND expressed in >= 20% samples

### Step B: STAR counts for limma
cd scripts
python download_gdc_counts.py
--subtype-csv ../01_Figure1_clustering/data/subtype_assignment_balanced.csv

text

First run takes 10-30 min; subsequent runs are fast because TSVs are
cached under `data/tcga_kirc/_gdc_cache/`.

## Inputs required (place them yourself)

### For `tcga_kirc_prep.py`

| File | Place at | Source |
|---|---|---|
| `TCGA-KIRC_gene_annotation.csv` | `data/raw/` | GDC / Xena download |
| `mrna/TCGA-KIRC_tpm_mrna.csv` | `data/raw/` | GDC / Xena download |
| `TCGA-KIRC_clinical_survival.csv` | `data/raw/` | GDC clinical export |
| `TCGA-KIRC_clinical_info.csv` | `data/raw/` | GDC clinical export |

### For `download_gdc_counts.py`

| File | Place at | Source |
|---|---|---|
| `subtype_assignment_balanced.csv` | `data/processed/` | `01_Figure1_clustering` |

## Manual copy

```bash
# Required before download_gdc_counts.py
cp ../01_Figure1_clustering/data/subtype_assignment_balanced.csv \
   data/processed/subtype_assignment_balanced.csv
Raw TCGA downloads are user-supplied; place the four files under
data/raw/.

Outputs
From tcga_kirc_prep.py
File	Consumed by
data/tcga_kirc/KIRC_expr_log2_tpm.csv	01, 02, 03, 04, 05, 06, 07, 08
data/tcga_kirc/KIRC_expr_tpm.csv	04 (xCell, CIBERSORT)
data/tcga_kirc/KIRC_clinical_cleaned.csv	01, 05, 06
data/tcga_kirc/KIRC_data_summary.json	informational
From download_gdc_counts.py
File	Consumed by
data/tcga_kirc/KIRC_counts_for_limma.csv	02, 04
data/tcga_kirc/sample_info_for_limma.csv	02
data/tcga_kirc/_gdc_cache/*.tsv	intermediate cache
Environment
Python >= 3.10: pandas, numpy, requests.

Notes
Both scripts auto-detect the project root by walking up the directory
tree until they find a folder named data/. This makes them robust to
the project folder being renamed.

tcga_kirc_prep.py accepts --data-dir and --output-dir, plus the
environment variables TCGA_KIRC_RAW and NCC_DATA_ROOT.

download_gdc_counts.py requires internet access to
api.gdc.cancer.gov.

The two entry scripts are independent.