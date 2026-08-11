"""Seaborn pairplot (scatterplot matrix) of the three tau variants: unweighted,
GLS-weighted, and ridge-regularized GLS-weighted."""
import sys
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use("Agg")
from plot_utils import annotate_diag_stats

IN_CSV = sys.argv[1]
OUT_PNG = sys.argv[2]

COLS = ["tau_unweighted", "tau_weighted_gls", "tau_weighted_gls_ridge"]

df = pd.read_csv(IN_CSV)[COLS].dropna()

g = sns.pairplot(df, plot_kws={"s": 5, "alpha": 0.15})
for i, yi in enumerate(COLS):
    for j, xj in enumerate(COLS):
        if i == j:
            continue
        ax = g.axes[i, j]
        ax.plot([0, 1], [0, 1], "--", color="gray", lw=1)
        annotate_diag_stats(ax, df[xj], df[yi], loc=(0.05, 0.95), fontsize=6)
g.figure.suptitle(f"tau variants, n={len(df)} genes", y=1.02)
g.savefig(OUT_PNG, dpi=300)
