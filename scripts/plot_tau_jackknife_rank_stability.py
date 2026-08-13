"""Jackknife counterpart to plot_tau_bootstrap_rank_stability.py: same layout, but
each point is one leave-one-unit-out rep (tissue, or whole t=10/t=20 cluster block)
instead of one bootstrap resample. Points above y=x mean ridge preserved gene ranking
better than unweighted when that unit was dropped."""
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plot_utils import annotate_diag_stats

IN_UNWEIGHTED_TISSUE = sys.argv[1]
IN_RIDGE_TISSUE = sys.argv[2]
IN_UNWEIGHTED_CLUSTER10 = sys.argv[3]
IN_RIDGE_CLUSTER10 = sys.argv[4]
IN_UNWEIGHTED_CLUSTER20 = sys.argv[5]
IN_RIDGE_CLUSTER20 = sys.argv[6]
OUT_PNG = sys.argv[7]
SUBSET_LABEL = sys.argv[8] if len(sys.argv) > 8 else None  # e.g. "highspec"


def load_pair(unweighted_csv, ridge_csv, mode):
    ucol, rcol = f"spearman_rho_unweighted_{mode}", f"spearman_rho_ridge_{mode}"
    u = pd.read_csv(unweighted_csv)[["rep", ucol]]
    w = pd.read_csv(ridge_csv)[["rep", rcol]]
    return u.merge(w, on="rep", how="inner").dropna(), ucol, rcol


fig, axes = plt.subplots(1, 3, figsize=(19, 5.5))

for ax, unweighted_csv, ridge_csv, mode, title in [
    (axes[0], IN_UNWEIGHTED_TISSUE, IN_RIDGE_TISSUE, "tissue", "leave-one-tissue-out"),
    (axes[1], IN_UNWEIGHTED_CLUSTER10, IN_RIDGE_CLUSTER10, "cluster10", "leave-one-block-out (t=10)"),
    (axes[2], IN_UNWEIGHTED_CLUSTER20, IN_RIDGE_CLUSTER20, "cluster20", "leave-one-block-out (t=20)"),
]:
    df, ucol, rcol = load_pair(unweighted_csv, ridge_csv, mode)
    ax.scatter(df[ucol], df[rcol], s=14, alpha=0.6, color="tab:orange")
    lo = min(df[ucol].min(), df[rcol].min())
    ax.plot([lo, 1], [lo, 1], "--", color="gray", lw=1)
    ax.set_xlabel("spearman rho: tau_unweighted (leave-out vs canonical)")
    ax.set_ylabel("spearman rho: tau_weighted_gls_ridge (leave-out vs canonical)")
    ax.set_title(f"{title} (n={len(df)} units)", fontsize=10)
    annotate_diag_stats(ax, df[ucol], df[rcol])
    print(f"{mode}: n={len(df)} median_unweighted={df[ucol].median():.4f} "
          f"median_ridge={df[rcol].median():.4f}", file=sys.stderr)

title = ("Rank stability of tau under jackknifing: unweighted vs ridge-weighted "
         "(Spearman rho vs canonical tau ranking, per leave-out unit)")
if SUBSET_LABEL:
    title += f" -- gene subset: {SUBSET_LABEL}"
fig.suptitle(title)
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)
