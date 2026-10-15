"""PERMANOVA sensitivity check for Figure 3 Panel B.

Panel B runs PERMANOVA on the significant DEGs, which were themselves
selected by the S1/S2 contrast. To rule out a selection artefact, this
script reruns the same PERMANOVA on all genes that survived the edgeR
filterByExpr step (i.e. the full set of rows in deg_limma_voom.csv).

Saves results/tables/permanova_allgenes_sensitivity_v1.csv.
"""

import os, sys, io
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from sklearn.preprocessing import StandardScaler
from sklearn.utils import resample

PROJ_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_ROOT = os.path.join(PROJ_ROOT, 'data')
RAW_DIR = os.path.join(DATA_ROOT, 'raw')
PROC_DIR = os.path.join(DATA_ROOT, 'processed')

if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8',
                                  errors='replace')


def permanova_anderson(X, groups, n_perm=999, seed=42):
    # Anderson 2001 pseudo-F on squared Euclidean distances,
    # permutation test with n_perm shuffles.
    D_sq = squareform(pdist(X, 'sqeuclidean'))
    N = X.shape[0]
    n1 = int(groups.sum())
    n2 = N - n1

    SS_total = D_sq.sum() / N
    SS_within = (D_sq[groups][:, groups].sum() / n1
                 + D_sq[~groups][:, ~groups].sum() / n2)
    SS_between = SS_total - SS_within

    f_obs = (SS_between / (2 - 1)) / (SS_within / (N - 2))
    R2 = SS_between / SS_total

    rng = np.random.RandomState(seed)
    idx = np.arange(N)
    fp = np.zeros(n_perm)
    for i in range(n_perm):
        pi = resample(idx, replace=False, n_samples=N, random_state=rng)
        p1 = pi[:n1]
        p2 = pi[n1:]
        W = (D_sq[p1][:, p1].sum() / n1
             + D_sq[p2][:, p2].sum() / n2)
        B = D_sq.sum() / N - W
        fp[i] = (B / (2 - 1)) / (W / (N - 2))

    p_perm = (np.sum(fp >= f_obs) + 1) / (n_perm + 1)
    return f_obs, R2, p_perm


deg_limma = pd.read_csv(os.path.join(PROC_DIR, 'deg_limma_voom.csv')
                        ).dropna(subset=['logFC'])
expr = pd.read_csv(os.path.join(RAW_DIR, 'KIRC_expr_log2_tpm.csv'),
                   index_col=0)
sub = pd.read_csv(os.path.join(PROC_DIR, 'subtype_assignment_balanced.csv'))
counts = pd.read_csv(os.path.join(RAW_DIR, 'KIRC_counts_for_limma.csv'),
                     index_col=0)

sub_map = dict(zip(sub['Patient'], sub['Subtype']))
s1_pats = [p for p, v in sub_map.items() if v == 'Subtype_1']
s2_pats = [p for p, v in sub_map.items() if v == 'Subtype_2']
common_samples = [c for c in expr.columns if c in s1_pats or c in s2_pats]
print(f"Samples: {len(common_samples)} ({len(s1_pats)} S1 + "
      f"{len(s2_pats)} S2)")

counts = counts.reindex(columns=common_samples)
cpm = counts.div(counts.sum(axis=0), axis=1) * 1e6
expr_ensg = np.log2(cpm + 1)

groups = np.array([sub_map.get(s) == 'Subtype_1' for s in common_samples])

# Run 1: DEG space (should reproduce the manuscript value)
deg_ensg_all = deg_limma[
    deg_limma['category'].isin(['UP_in_S1', 'DOWN_in_S1'])
]['gene'].tolist()
deg_genes = [g for g in deg_ensg_all if g in expr_ensg.index]
X_deg = StandardScaler().fit_transform(
    expr_ensg.loc[deg_genes, common_samples].T.values)
f_deg, r2_deg, p_deg = permanova_anderson(X_deg, groups)
print(f"\nDEG space, n={len(deg_genes)}:  "
      f"pseudo-F={f_deg:.2f}  R2={r2_deg:.3f}  p={p_deg:.3f}")

# Run 2: all genes that survived filterByExpr
tested_genes = [g for g in deg_limma['gene'].tolist()
                if g in expr_ensg.index]
X_all = StandardScaler().fit_transform(
    expr_ensg.loc[tested_genes, common_samples].T.values)
f_all, r2_all, p_all = permanova_anderson(X_all, groups)
print(f"All tested, n={len(tested_genes)}:  "
      f"pseudo-F={f_all:.2f}  R2={r2_all:.3f}  p={p_all:.3f}")

res = pd.DataFrame([
    {'space': 'DEG (FDR<0.05, |log2FC|>0.5)',
     'n_genes': len(deg_genes),
     'pseudo_F': f_deg, 'R2': r2_deg, 'P_perm': p_deg},
    {'space': 'All tested genes (filterByExpr)',
     'n_genes': len(tested_genes),
     'pseudo_F': f_all, 'R2': r2_all, 'P_perm': p_all},
])
out_csv = os.path.join(PROJ_ROOT, 'results', 'tables',
                       'permanova_allgenes_sensitivity_v1.csv')
os.makedirs(os.path.dirname(out_csv), exist_ok=True)
res.to_csv(out_csv, index=False)
print(f"\nSaved: {out_csv}")