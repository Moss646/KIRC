"""
xCell Immune Infiltration Analysis for TCGA KIRC
=================================================
Uses pre-computed xCell scores (R xCell v1.1.0, spillover corrected).
Generates a 2x2 panel:
  A. Differential xCell scores (S1 minus S2)
  B. Key immune cell boxplots by functional group (activation / suppression / APC)
  C. xCell composite ImmuneScore vs StromaScore
  D. Heatmap of top-variant cell types
"""

import os
import sys
import io
import warnings
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

warnings.filterwarnings('ignore')

# ============================================================
# Paths
# ============================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJ_ROOT = os.path.dirname(SCRIPT_DIR)
DATA_ROOT = os.environ.get('MEL_DATA_ROOT', PROJ_ROOT)
TCGA_DIR = os.path.join(DATA_ROOT, 'data', 'tcga_kirc')
CONS_DIR = os.path.join(DATA_ROOT, 'data', 'processed')
XCELL_DIR = os.path.join(DATA_ROOT, 'data', 'xcell')
OUT_DIR = os.path.join(PROJ_ROOT, 'results', 'figures')
os.makedirs(OUT_DIR, exist_ok=True)

if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

# ============================================================
# Aesthetics
# ============================================================
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 14
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['ytick.major.width'] = 2

# ============================================================
# 1. Load data
# ============================================================
print('[xCell] Loading expression + subtype + xCell scores...')
EXPR = pd.read_csv(os.path.join(TCGA_DIR, 'KIRC_expr_log2_tpm.csv'), index_col=0)
EXPR_TPM = pd.read_csv(os.path.join(TCGA_DIR, 'KIRC_expr_tpm.csv'), index_col=0)
print(f'  Expression matrix: {EXPR_TPM.shape[0]} genes x {EXPR_TPM.shape[1]} samples')

sub = pd.read_csv(os.path.join(CONS_DIR, 'subtype_assignment_balanced.csv'))
sub_map = dict(zip(sub['Patient'], sub['Subtype']))
s1_samples = [p for p, v in sub_map.items() if v == 'Subtype_1']
s2_samples = [p for p, v in sub_map.items() if v == 'Subtype_2']
print(f'  S1: {len(s1_samples)}, S2: {len(s2_samples)}')

# ============================================================
# 2. xCell scores (R pre-computed, spillover corrected)
# ============================================================
xcell_scores = pd.read_csv(os.path.join(XCELL_DIR, 'xcell_scores_raw.csv'), index_col=0)
print(f'  xCell raw scores: {xcell_scores.shape[0]} cell types x {xcell_scores.shape[1]} samples')

xcell_z = pd.read_csv(os.path.join(XCELL_DIR, 'xcell_scores_zscore.csv'), index_col=0).T
print(f'  xCell z-scores: {xcell_z.shape} (samples x cell types)')

diff_df = pd.read_csv(os.path.join(XCELL_DIR, 'xcell_subtype_diff.csv'))
print(f'  Diff stats: {len(diff_df)} cell types, {(diff_df["fdr"] < 0.05).sum()} significant (FDR<0.05)')

# ============================================================
# 3. Subtype-stratified xCell
# ============================================================
common = [s for s in xcell_z.index if s in sub_map]
xcell_filt = xcell_z.loc[common].astype(float)
s1_idx = [s for s in common if sub_map[s] == 'Subtype_1']
s2_idx = [s for s in common if sub_map[s] == 'Subtype_2']

# ============================================================
# 4. Plot
# ============================================================
print('[xCell] Generating figure...')

# Panel A: top 30 significant CELL-TYPE scores only — composite scores
# (ImmuneScore / StromaScore / MicroenvironmentScore) are excluded here and
# shown separately in Panel C.
COMPOSITES = {'ImmuneScore', 'StromaScore', 'MicroenvironmentScore'}
_sig = diff_df[(diff_df['fdr'] < 0.05) & (~diff_df['cell_type'].isin(COMPOSITES))]
plot_cts = _sig.head(30)['cell_type'].tolist()
if len(plot_cts) < 15:
    plot_cts = _sig.head(20)['cell_type'].tolist()

fig, axes = plt.subplots(2, 2, figsize=(24, 22))
plt.subplots_adjust(hspace=0.35)

# ----- A: Differential xCell scores (S1 minus S2) -----
ax = axes[0, 0]
plot_data = diff_df[diff_df['cell_type'].isin(plot_cts)].copy()
plot_data['diff_s1s2'] = -plot_data['mean_zscore_diff']   # S1 minus S2 (column mean_zscore_diff = S2 minus S1)
plot_data = plot_data.sort_values('diff_s1s2', ascending=True)
colors = ['#E74C3C' if v > 0 else '#3498DB' for v in plot_data['diff_s1s2']]
ax.barh(range(len(plot_data)), plot_data['diff_s1s2'], color=colors, edgecolor='white', linewidth=0.5)
ax.set_yticks(range(len(plot_data)))
ax.set_yticklabels(plot_data['cell_type'], fontsize=12)

# significance stars at the tip of each bar (FDR from diff_df)
_xmin, _xmax = plot_data['diff_s1s2'].min(), plot_data['diff_s1s2'].max()
_xpad = (_xmax - _xmin) * 0.04
for _i, (_idx, _row) in enumerate(plot_data.iterrows()):
    _fdr = _row['fdr']
    _stars = '***' if _fdr < 0.001 else '**' if _fdr < 0.01 else '*' if _fdr < 0.05 else 'ns'
    _tip = _row['diff_s1s2']
    _sx = _tip + _xpad if _tip >= 0 else _tip - _xpad
    _ha = 'left' if _tip >= 0 else 'right'
    ax.text(_sx, _i, _stars, va='center', ha=_ha, fontsize=12,
            fontweight='bold', color='black')
ax.set_xlim(_xmin - _xpad * 2, _xmax + _xpad * 2.5)

ax.set_xlabel('Difference in xCell enrichment scores between S1 and S2 tumors (S1 minus S2)', fontsize=20)
ax.set_title('A', fontsize=26, fontweight='bold', loc='left', pad=8)
ax.axvline(0, color='black', linewidth=0.8)
legend_a = [Patch(facecolor='#E74C3C', alpha=0.85, label='Higher in S1'),
            Patch(facecolor='#3498DB', alpha=0.85, label='Higher in S2')]
ax.legend(handles=legend_a, fontsize=16, loc='upper left', framealpha=0.9)

# ----- B: Key immune cell boxplots (horizontal, paired, Fig4-Panel-D style) -----
ax = axes[0, 1]
# Display order (top -> bottom in the figure), grouped by biological function:
#   effector / cytotoxic lymphocytes -> antigen-presenting cells -> suppressive cells
panel_b_display = [
    ('Th1 cells', 'Th1 cells'),
    ('NKT', 'NKT cells'),
    ('CD8+ T-cells', 'CD8+ T cells'),
    ('B-cells', 'B cells'),
    ('DC', 'DC'),
    ('Tregs', 'Tregs'),
    ('Macrophages', 'Macrophages'),
]
# matplotlib draws y=0 at the bottom, so reverse for plotting and keep
# boxes / stars / tick labels all driven by the SAME list order (y = i).
panel_b_plot = list(reversed(panel_b_display))
b_cts = [d for d, _lab in panel_b_plot]
b_labels = [lab for _d, lab in panel_b_plot]

bp_data, bp_positions, bp_colors = [], [], []
for i, ct in enumerate(b_cts):
    s1_vals = xcell_filt.loc[s1_idx, ct].dropna().values
    s2_vals = xcell_filt.loc[s2_idx, ct].dropna().values
    bp_data.append(s1_vals); bp_positions.append(i - 0.15); bp_colors.append('#E74C3C')
    bp_data.append(s2_vals); bp_positions.append(i + 0.15); bp_colors.append('#3498DB')

bp = ax.boxplot(bp_data, positions=bp_positions, vert=False, patch_artist=True,
                widths=0.22, showfliers=False,
                medianprops=dict(color='black', linewidth=1.8))
for patch, color in zip(bp['boxes'], bp_colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.75)
    patch.set_edgecolor('white')

# group separators (bottom->top: suppression | antigen presentation | activation)
for ysep in [1.5, 3.5]:
    ax.axhline(ysep, color='gray', linestyle='--', linewidth=1, alpha=0.5)

# significance stars just right of the boxes (FDR from diff_df)
all_vals = np.concatenate([d for d in bp_data if len(d) > 0])

def _upper_whisker(v):
    q1, q3 = np.percentile(v, [25, 75])
    return q3 + 1.5 * (q3 - q1)

star_x = max(_upper_whisker(d) for d in bp_data if len(d) > 0) + 0.5
for i, ct in enumerate(b_cts):
    fdr = diff_df.loc[diff_df['cell_type'] == ct, 'fdr']
    fdr = fdr.values[0] if len(fdr) else 1.0
    stars = ('***' if fdr < 0.001 else '**' if fdr < 0.01 else '*' if fdr < 0.05 else 'ns')
    ax.text(star_x, i, stars, va='center', ha='left', fontsize=18,
            fontweight='bold', color='black')

ax.set_yticks(range(len(b_cts)))
ax.set_yticklabels(b_labels, fontsize=16)
ax.set_xlabel('Row-wise Z-score of xCell enrichment scores', fontsize=20)
ax.set_title('B', fontsize=26, fontweight='bold', loc='left', pad=8)
ax.set_ylim(-0.8, len(b_cts) - 0.2)
ax.set_xlim(-3, 4)
ax.set_xticks(np.arange(-2, 5, 2))
legend_elements = [Patch(facecolor='#E74C3C', alpha=0.75, label='S1'),
                   Patch(facecolor='#3498DB', alpha=0.75, label='S2')]
ax.legend(handles=legend_elements, fontsize=14, loc='upper right')

# ----- C: xCell composite ImmuneScore vs StromaScore -----
ax = axes[1, 0]
# Use xCell's native composite scores (one value per sample) instead of a
# custom mean over hand-picked cell-type lists.
comp = xcell_scores.loc[['StromaScore', 'ImmuneScore']]  # rows = score, cols = samples
s1c = [c for c in s1_samples if c in comp.columns]
s2c = [c for c in s2_samples if c in comp.columns]
ax.scatter(comp.loc['StromaScore', s1c], comp.loc['ImmuneScore', s1c],
           c='#E74C3C', label='S1', alpha=0.5, s=60, marker='o')
ax.scatter(comp.loc['StromaScore', s2c], comp.loc['ImmuneScore', s2c],
           c='#3498DB', label='S2', alpha=0.5, s=60, marker='^')
ax.set_xlabel('xCell StromaScore', fontsize=20)
ax.set_ylabel('xCell ImmuneScore', fontsize=20)
ax.set_title('C', fontsize=26, fontweight='bold', loc='left', pad=8)
ax.legend(fontsize=14, loc='upper right', markerscale=2)
ax.tick_params(labelsize=18)
# keep variable names referenced by the summary print below
immune_score_s1 = comp.loc['ImmuneScore', s1c]
immune_score_s2 = comp.loc['ImmuneScore', s2c]

# ----- D: Heatmap -----
ax = axes[1, 1]
heatmap_data = pd.DataFrame({
    'S1': xcell_filt.loc[s1_idx].mean(),
    'S2': xcell_filt.loc[s2_idx].mean(),
}).T
variances = xcell_filt.var()
top_var = variances.nlargest(25).index
heatmap_data = heatmap_data[top_var]

im = ax.imshow(heatmap_data.values, aspect='auto', cmap='RdBu_r', vmin=-1, vmax=1)
ax.set_xticks(range(len(top_var)))
ax.set_xticklabels(top_var, fontsize=18, rotation=90, ha='center')
ax.set_yticks(range(2))
ax.set_yticklabels(['S1', 'S2'], fontsize=18)
ax.set_title('D', fontsize=26, fontweight='bold', loc='left', pad=8)
cbar = plt.colorbar(im, ax=ax, shrink=0.8)
cbar.set_label('Mean z-score', fontsize=20)
cbar.ax.tick_params(labelsize=18)

plt.tight_layout(pad=3)
out_png = os.path.join(OUT_DIR, 'Supplementary_xCell_analysis_v18_compfilter.png')
out_svg = os.path.join(OUT_DIR, 'Supplementary_xCell_analysis_v18_compfilter.svg')
fig.savefig(out_png, dpi=200, bbox_inches='tight', facecolor='white')
fig.savefig(out_svg, format='svg', bbox_inches='tight', facecolor='white')
plt.close()

print(f'[xCell] Saved: {out_png}')
print(f'[xCell] Saved: {out_svg}')
print(f'[xCell] S1 mean immune score: {immune_score_s1.mean():.3f}')
print(f'[xCell] S2 mean immune score: {immune_score_s2.mean():.3f}')
print()
print('Top 10 differential cell types (S2 vs S1):')
for _, row in diff_df.head(10).iterrows():
    stars = '***' if row['fdr'] < 0.001 else '**' if row['fdr'] < 0.01 else '*' if row['fdr'] < 0.05 else 'ns'
    direction = '↑S2' if row['mean_zscore_diff'] > 0 else '↓S2'
    print(f'  {row["cell_type"]:<30s} mean_zscore_diff={row["mean_zscore_diff"]:+7.3f}  FDR={row["fdr"]:.2e} {stars} {direction}')
