"""3x3 extension of plot_tau_jackknife_rank_stability.py: columns are the three
leave-out schemes (tissue / cluster10 / cluster20); rows compare tau_unweighted
against ridge-weighted, cluster-size-downweighted (k=10), and cluster-size-downweighted
(k=20) tau (see compute_tau_clustersize.py). Points above y=x mean the alternative
weighting preserved gene ranking better than unweighted when that unit was dropped."""
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
IN_CLUSTERSIZE10_TISSUE = sys.argv[7]
IN_CLUSTERSIZE10_CLUSTER10 = sys.argv[8]
IN_CLUSTERSIZE10_CLUSTER20 = sys.argv[9]
IN_CLUSTERSIZE20_TISSUE = sys.argv[10]
IN_CLUSTERSIZE20_CLUSTER10 = sys.argv[11]
IN_CLUSTERSIZE20_CLUSTER20 = sys.argv[12]
OUT_PNG = sys.argv[13]

MODE_TITLE = {
    "tissue": "leave-one-tissue-out",
    "cluster10": "leave-one-block-out (t=10)",
    "cluster20": "leave-one-block-out (t=20)",
}


def load_pair(unweighted_csv, y_csv, ycol, mode):
    ucol = f"spearman_rho_unweighted_{mode}"
    u = pd.read_csv(unweighted_csv)[["rep", ucol]]
    y = pd.read_csv(y_csv)[["rep", ycol]]
    return u.merge(y, on="rep", how="inner").dropna(), ucol


ROWS = [
    ("tau_weighted_gls_ridge", "spearman_rho_ridge_{mode}", "tab:orange", [
        (IN_UNWEIGHTED_TISSUE, IN_RIDGE_TISSUE, "tissue"),
        (IN_UNWEIGHTED_CLUSTER10, IN_RIDGE_CLUSTER10, "cluster10"),
        (IN_UNWEIGHTED_CLUSTER20, IN_RIDGE_CLUSTER20, "cluster20"),
    ]),
    ("tau_clustersize (k=10)", "spearman_rho_clustersize_{mode}", "tab:green", [
        (IN_UNWEIGHTED_TISSUE, IN_CLUSTERSIZE10_TISSUE, "tissue"),
        (IN_UNWEIGHTED_CLUSTER10, IN_CLUSTERSIZE10_CLUSTER10, "cluster10"),
        (IN_UNWEIGHTED_CLUSTER20, IN_CLUSTERSIZE10_CLUSTER20, "cluster20"),
    ]),
    ("tau_clustersize (k=20)", "spearman_rho_clustersize_{mode}", "tab:purple", [
        (IN_UNWEIGHTED_TISSUE, IN_CLUSTERSIZE20_TISSUE, "tissue"),
        (IN_UNWEIGHTED_CLUSTER10, IN_CLUSTERSIZE20_CLUSTER10, "cluster10"),
        (IN_UNWEIGHTED_CLUSTER20, IN_CLUSTERSIZE20_CLUSTER20, "cluster20"),
    ]),
]

fig, axes = plt.subplots(3, 3, figsize=(19, 16.5))

for row_i, (y_label, ycol_tmpl, color, panels) in enumerate(ROWS):
    for col_i, (unweighted_csv, y_csv, mode) in enumerate(panels):
        ax = axes[row_i, col_i]
        ycol = ycol_tmpl.format(mode=mode)
        df, ucol = load_pair(unweighted_csv, y_csv, ycol, mode)
        ax.scatter(df[ucol], df[ycol], s=14, alpha=0.6, color=color)
        lo = min(df[ucol].min(), df[ycol].min())
        ax.plot([lo, 1], [lo, 1], "--", color="gray", lw=1)
        ax.set_xlabel("spearman rho: tau_unweighted (leave-out vs canonical)")
        ax.set_ylabel(f"spearman rho: {y_label} (leave-out vs canonical)")
        title = MODE_TITLE[mode] if row_i == 0 else ""
        ax.set_title(f"{title} (n={len(df)} units)" if title else f"n={len(df)} units", fontsize=10)
        annotate_diag_stats(ax, df[ucol], df[ycol])
        print(f"{y_label}/{mode}: n={len(df)} median_unweighted={df[ucol].median():.4f} "
              f"median_y={df[ycol].median():.4f}", file=sys.stderr)

fig.suptitle("Rank stability of tau under jackknifing: unweighted vs ridge-weighted / "
             "cluster-size-downweighted (Spearman rho vs canonical tau ranking, per leave-out unit)")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)
