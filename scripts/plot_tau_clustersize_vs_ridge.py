"""Pairplot comparing tau_unweighted, tau_weighted_gls_ridge (default GLS/ridge
pipeline) against tau_clustersize at k=10 and k=20 (cluster-size downweighting, see
compute_tau_clustersize.py) -- direct comparison of the two tissue-weighting
approaches."""
import sys
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use("Agg")
from plot_utils import annotate_diag_stats

IN_RIDGE_TAU = sys.argv[1]
IN_CLUSTERSIZE10_TAU = sys.argv[2]
IN_CLUSTERSIZE20_TAU = sys.argv[3]
OUT_PNG = sys.argv[4]

ridge = pd.read_csv(IN_RIDGE_TAU)[["gene", "tau_unweighted", "tau_weighted_gls_ridge"]]
c10 = pd.read_csv(IN_CLUSTERSIZE10_TAU)[["gene", "tau_clustersize"]].rename(
    columns={"tau_clustersize": "tau_clustersize_k10"})
c20 = pd.read_csv(IN_CLUSTERSIZE20_TAU)[["gene", "tau_clustersize"]].rename(
    columns={"tau_clustersize": "tau_clustersize_k20"})

df = ridge.merge(c10, on="gene").merge(c20, on="gene")
COLS = ["tau_unweighted", "tau_weighted_gls_ridge", "tau_clustersize_k10", "tau_clustersize_k20"]
df = df[COLS].dropna()

g = sns.pairplot(df, plot_kws={"s": 5, "alpha": 0.15})
for i, yi in enumerate(COLS):
    for j, xj in enumerate(COLS):
        if i == j:
            continue
        ax = g.axes[i, j]
        ax.plot([0, 1], [0, 1], "--", color="gray", lw=1)
        annotate_diag_stats(ax, df[xj], df[yi], loc=(0.05, 0.95), fontsize=6)
g.figure.suptitle(f"GLS/ridge vs cluster-size-downweighted tau, n={len(df)} genes", y=1.02)
g.savefig(OUT_PNG, dpi=300)
