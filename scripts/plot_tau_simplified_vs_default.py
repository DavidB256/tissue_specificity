"""3-panel comparison of tau on default (68) vs simplified (fine-grained sub-tissues
collapsed into parent, see build_gct_simplified.py) tissues: (left) tau_unweighted,
(center) tau_weighted_gls_ridge with the default (unsquared) corr mat, (right)
tau_weighted_gls_ridge with the elementwise-squared corr mat (see
compute_weighted_tau_squared.py)."""
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plot_utils import annotate_diag_stats

IN_DEFAULT_TAU = sys.argv[1]
IN_DEFAULT_TAU_SQUARED = sys.argv[2]
IN_SIMPLIFIED_TAU = sys.argv[3]
IN_SIMPLIFIED_TAU_SQUARED = sys.argv[4]
OUT_PNG = sys.argv[5]

d = pd.read_csv(IN_DEFAULT_TAU)[["gene", "tau_unweighted", "tau_weighted_gls_ridge"]]
d_sq = pd.read_csv(IN_DEFAULT_TAU_SQUARED)[["gene", "tau_weighted_gls_ridge"]].rename(
    columns={"tau_weighted_gls_ridge": "tau_weighted_gls_ridge_sq"})
s = pd.read_csv(IN_SIMPLIFIED_TAU)[["gene", "tau_unweighted", "tau_weighted_gls_ridge"]].rename(
    columns={"tau_unweighted": "tau_unweighted_simp", "tau_weighted_gls_ridge": "tau_weighted_gls_ridge_simp"})
s_sq = pd.read_csv(IN_SIMPLIFIED_TAU_SQUARED)[["gene", "tau_weighted_gls_ridge"]].rename(
    columns={"tau_weighted_gls_ridge": "tau_weighted_gls_ridge_simp_sq"})

df = d.merge(d_sq, on="gene").merge(s, on="gene").merge(s_sq, on="gene").dropna()

PANELS = [
    ("tau_unweighted", "tau_unweighted_simp", "tau_unweighted: default vs simplified tissues"),
    ("tau_weighted_gls_ridge", "tau_weighted_gls_ridge_simp",
     "tau_weighted_gls_ridge (default R): default vs simplified tissues"),
    ("tau_weighted_gls_ridge_sq", "tau_weighted_gls_ridge_simp_sq",
     "tau_weighted_gls_ridge (squared R): default vs simplified tissues"),
]

AXIS_LO, AXIS_HI = 0.2, 1.0

fig, axes = plt.subplots(1, 3, figsize=(18, 6))
for ax, (xcol, ycol, title) in zip(axes, PANELS):
    ax.scatter(df[xcol], df[ycol], s=4, alpha=0.1, color="tab:blue")
    ax.plot([AXIS_LO, AXIS_HI], [AXIS_LO, AXIS_HI], "--", color="gray", lw=1)
    ax.set_xlim(AXIS_LO, AXIS_HI)
    ax.set_ylim(AXIS_LO, AXIS_HI)
    ax.set_xlabel(f"{xcol} (default, 68 tissues)")
    ax.set_ylabel(f"{ycol} (simplified tissues)")
    ax.set_title(title, fontsize=10)
    annotate_diag_stats(ax, df[xcol], df[ycol])

fig.suptitle(f"tau: default vs simplified tissues, n={len(df)} genes")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)
