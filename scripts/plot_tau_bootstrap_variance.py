"""Visualize per-gene bootstrap variance of tau_unweighted vs tau_weighted_gls_ridge,
under (i) tissue resampling and (ii) correlation-block resampling. Log-log scatter
since variances span orders of magnitude across genes."""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plot_utils import annotate_diag_stats

IN_UNWEIGHTED_TISSUE = sys.argv[1]
IN_RIDGE_TISSUE = sys.argv[2]
IN_UNWEIGHTED_CLUSTER = sys.argv[3]
IN_RIDGE_CLUSTER = sys.argv[4]
OUT_PNG = sys.argv[5]


def load_pair(unweighted_csv, ridge_csv, mode):
    ucol, rcol = f"tau_unweighted_{mode}_var", f"tau_ridge_{mode}_var"
    u = pd.read_csv(unweighted_csv)[["gene", ucol]]
    w = pd.read_csv(ridge_csv)[["gene", rcol]]
    df = u.merge(w, on="gene", how="inner").replace([np.inf, -np.inf], np.nan).dropna()
    df = df[(df[ucol] > 0) & (df[rcol] > 0)]
    return df, ucol, rcol


fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))

for ax, unweighted_csv, ridge_csv, mode, title in [
    (axes[0], IN_UNWEIGHTED_TISSUE, IN_RIDGE_TISSUE, "tissue", "tissue resampling"),
    (axes[1], IN_UNWEIGHTED_CLUSTER, IN_RIDGE_CLUSTER, "cluster", "correlation-block resampling"),
]:
    df, ucol, rcol = load_pair(unweighted_csv, ridge_csv, mode)
    x, y = df[ucol], df[rcol]
    ax.scatter(x, y, s=4, alpha=0.15, color="tab:blue")
    lo, hi = min(x.min(), y.min()), max(x.max(), y.max())
    ax.plot([lo, hi], [lo, hi], "--", color="gray", lw=1)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Var(tau_unweighted) across resamples")
    ax.set_ylabel("Var(tau_weighted_gls_ridge) across resamples")
    ax.set_title(f"{title} (n={len(df)} genes)", fontsize=10)
    annotate_diag_stats(ax, x, y)
    print(f"{mode}: n={len(df)}", file=sys.stderr)

fig.suptitle("Bootstrap variance of tau: unweighted vs ridge-weighted")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)
