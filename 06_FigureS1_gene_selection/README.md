# 06 — Gene-set collection and 68-gene panel selection

Builds the lipid-metabolism gene universe from MSigDB (Hallmark, KEGG,
Reactome, GO:BP), hierarchically de-duplicates the GO:BP sets by information
content and ancestor coverage, screens every gene by univariate Cox, applies
strict multi-criteria selection, and assembles the final 68-gene balanced
panel. Produces the Figure 1 overview plot.

## Pipeline
cd scripts
python run_upstream.py # steps 1-5 below, then figure
python run_upstream.py --start 3 --end 5 --no-figure

text

| Step | Script | Purpose |
|---|---|---|
| 1 | `download_msigdb_lipid.py` | Download Hallmark / KEGG / Reactome / GO:BP GMTs, filter by lipid keywords, build the pool |
| 2 | `go_dedup_lipid.py` | Hierarchical GO:BP de-duplication by IC and ancestor coverage |
| 3 | `cox_screening.py` | Univariate Cox across all lipid genes |
| 4 | `strict_selection.py` | p < 0.001, \|log2HR\| > 0.30, SD > median, \|r\| < 0.85, layer assignment |
| 5 | `balance_genes_for_clustering.py` | Quota-based 68-gene panel (20 + 15 + 15 + 5 + 2 + 11) |
| 6 | `gene_selection_fig_v12_pubsize.py` | Figure 1 overview (volcano, top-15 HR, layer distribution) |

## Inputs required (place them yourself)

| File | Place at | Source |
|---|---|---|
| `KIRC_expr_log2_tpm.csv` | `data/tcga_kirc/` | `09_TCGA-KIRC_prep` |
| `KIRC_clinical_cleaned.csv` | `data/tcga_kirc/` | `09_TCGA-KIRC_prep` |
| MSigDB GMTs | (downloaded by step 1) | MSigDB release |
| `go-basic.obo` | (downloaded by step 2) | GO Consortium |

## Manual copy

```bash
cp ../09_TCGA-KIRC_prep/data/tcga_kirc/KIRC_expr_log2_tpm.csv     data/tcga_kirc/
cp ../09_TCGA-KIRC_prep/data/tcga_kirc/KIRC_clinical_cleaned.csv data/tcga_kirc/
Steps 1 and 2 download their own inputs from MSigDB and the GO Consortium.

Outputs (selected)
File	Consumed by
data/consensus_cluster/balanced_genes.txt	01_Figure1_clustering
data/msigdb_lipid_gene_pool_dedup.csv	01_Figure1_clustering (survival-blind)
data/msigdb_lipid_genesets_dedup.csv	02_Figure2_transcriptomic_CPTAC
data/gene_selection/selected_genes_strict.csv	04_Figure4_immune
results/Figure_1_gene_selection_v12_pubsize.png/.svg	manuscript Figure 1
Intermediate outputs (msigdb_lipid_genesets_summary.csv,
cox_screening/*, gene_selection/*, go_*.csv/json) are also produced
and used internally.

Environment
Python >= 3.10: pandas, numpy, scipy, lifelines, matplotlib,
requests.

Steps 1 and 2 require internet access. For offline runs, pre-place the
MSigDB GMTs in data/ and (optionally) the GO OBO file, and skip step 1
by running with --start 2.

Notes
run_upstream.py accepts --start, --end, --no-figure, and
--classifier. The optional --classifier argument points to a
reference subtype_classifier.json; if provided, the pipeline verifies
that the freshly selected 68-gene panel matches the reference set.

The final 68-gene panel is deterministic given the same input and
selection thresholds.