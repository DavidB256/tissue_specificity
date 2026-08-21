"""Cluster-weighted analog of plot_tau_simplified_vs_default.py: 2-panel comparison of
tau on default (68) vs a collapsed-tissue variant (LABEL, e.g. "onebrain") tissues, using
cluster-size-based weights instead of GLS (see compute_tau_cluster_weighted.py): (1)
tau_unweighted, (2) tau_weighted_cluster."""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plot_utils import annotate_diag_stats

IN_DEFAULT_TAU = sys.argv[1]
IN_COLLAPSED_TAU = sys.argv[2]
OUT_PNG = sys.argv[3]
LABEL = sys.argv[4] if len(sys.argv) > 4 else "onebrain"

d = pd.read_csv(IN_DEFAULT_TAU)[["gene", "tau_unweighted", "tau_weighted_cluster"]]
c = pd.read_csv(IN_COLLAPSED_TAU)[["gene", "tau_unweighted", "tau_weighted_cluster"]].rename(
    columns={"tau_unweighted": "tau_unweighted_simp", "tau_weighted_cluster": "tau_weighted_cluster_simp"})

df = d.merge(c, on="gene").dropna()

PANELS = [
    ("tau_unweighted", "tau_unweighted_simp", f"tau_unweighted: default vs {LABEL} tissues"),
    ("tau_weighted_cluster", "tau_weighted_cluster_simp",
     f"tau_weighted_cluster (1/cluster_size): default vs {LABEL} tissues"),
]

AXIS_LO, AXIS_HI = 0.2, 1.0

fig, axes = plt.subplots(1, 2, figsize=(12, 6))
for ax, (xcol, ycol, title) in zip(axes, PANELS):
    ax.scatter(df[xcol], df[ycol], s=4, alpha=0.1, color="tab:blue")
    ax.plot([AXIS_LO, AXIS_HI], [AXIS_LO, AXIS_HI], "--", color="gray", lw=1)
    ax.set_xlim(AXIS_LO, AXIS_HI)
    ax.set_ylim(AXIS_LO, AXIS_HI)
    ax.set_xlabel(f"{xcol} (default, 68 tissues)")
    ax.set_ylabel(f"{ycol} ({LABEL} tissues)")
    ax.set_title(title, fontsize=10)
    annotate_diag_stats(ax, df[xcol], df[ycol])

    dist = (df[xcol] - df[ycol]).abs() / np.sqrt(2)
    ax.legend(handles=[plt.Line2D([0], [0], color="none",
                                   label=f"dist from y=x: mean={dist.mean():.4f}, median={dist.median():.4f}")],
              loc="lower right", handlelength=0, handletextpad=0, fontsize=8)

fig.suptitle(f"cluster-weighted tau: default vs {LABEL} tissues, n={len(df)} genes")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)
