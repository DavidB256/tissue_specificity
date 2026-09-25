"""Variant of plot_tau_scenarios_scatter_thresholds.py that colors each gene by the
delete-one-cluster jackknife variance (see jackknife_cluster_variance.py) of whichever
tau estimate is on that panel's x-axis, instead of by threshold-crossing category.
Log color scale (variances span orders of magnitude); non-positive variances floored to
a small epsilon for display only."""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from plot_utils import annotate_diag_stats, display_label

IN_FREE_TAU = sys.argv[1]
IN_REF_TAU = sys.argv[2]
IN_JACKKNIFE_SUMMARY = sys.argv[3]
OUT_PNG = sys.argv[4]
LABEL = sys.argv[5] if len(sys.argv) > 5 else "reference"
UNIT = sys.argv[6] if len(sys.argv) > 6 else "cluster"

free = pd.read_csv(IN_FREE_TAU)[["gene", "tau_unweighted", "tau_free_weighted"]]
ref = pd.read_csv(IN_REF_TAU)[["gene", "tau_unweighted"]].rename(columns={"tau_unweighted": "tau_ref"})
jk = pd.read_csv(IN_JACKKNIFE_SUMMARY)[["gene", "tau_unweighted_jk_var", "tau_free_cluster_weighted_jk_var"]]
df = free.merge(ref, on="gene").merge(jk, on="gene").dropna()

PANELS = [
    ("tau_unweighted", "tau", "tau (full)", "tau_unweighted_jk_var"),
    ("tau_free_weighted", "tau_weighted", "tau_weighted (full, free/DE-optimized)", "tau_free_cluster_weighted_jk_var"),
]

fig, axes = plt.subplots(1, 2, figsize=(16, 7.5))
for ax, (xcol, title, xlabel, varcol) in zip(axes, PANELS):
    var = df[varcol].values
    positive = var[var > 0]
    vmin = np.quantile(positive, 0.1) if len(positive) else 1e-12
    vmax = np.quantile(positive, 0.95) if len(positive) else 1e-6
    var_plot = np.clip(var, vmin, vmax)
    sc = ax.scatter(df[xcol], df["tau_ref"], c=var_plot, s=6, alpha=0.6, cmap="viridis",
                     norm=LogNorm(vmin=vmin, vmax=vmax), rasterized=True)
    cb = fig.colorbar(sc, ax=ax)
    cb.set_label(f"jackknife variance ({varcol})")
    ax.plot([0, 1], [0, 1], "--", color="gray", lw=1, zorder=1)
    for thresh in (0.5, 0.8):
        ax.axvline(thresh, color="lightgray", lw=0.8, zorder=0)
        ax.axhline(thresh, color="lightgray", lw=0.8, zorder=0)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(f"tau_unweighted ({display_label(LABEL)})")
    ax.set_title(title, fontsize=11)
    annotate_diag_stats(ax, df[xcol], df["tau_ref"])

fig.suptitle(f"full vs {display_label(LABEL)}: tau and tau_weighted, colored by delete-one-{UNIT} "
             f"jackknife variance, n={len(df)} genes")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)

print(f"genes={len(df)}", file=sys.stderr)
print(f"wrote {OUT_PNG}", file=sys.stderr)
