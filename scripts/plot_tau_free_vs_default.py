"""free_weights analog of plot_tau_cluster_vs_default.py: 2-panel comparison of tau on
default (68 tissues) vs onebrain (56 tissues, unweighted reference used as the
optimization target, see optimize_free_weights.py): (1) tau_unweighted vs tau_unweighted,
(2) tau_free_weighted (default) vs tau_unweighted (onebrain)."""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plot_utils import annotate_diag_stats

IN_FREE_TAU = sys.argv[1]
IN_ONEBRAIN_TAU = sys.argv[2]
OUT_PNG = sys.argv[3]

d = pd.read_csv(IN_FREE_TAU)[["gene", "tau_unweighted", "tau_free_weighted"]]
c = pd.read_csv(IN_ONEBRAIN_TAU)[["gene", "tau_unweighted"]].rename(columns={"tau_unweighted": "tau_unweighted_ob"})

df = d.merge(c, on="gene").dropna()

PANELS = [
    ("tau_unweighted", "tau_unweighted_ob", "tau_unweighted"),
    ("tau_free_weighted", "tau_unweighted_ob", "tau_free_weighted"),
]

AXIS_LO, AXIS_HI = 0.2, 1.0

fig, axes = plt.subplots(1, 2, figsize=(12, 6))
for ax, (xcol, ycol, title) in zip(axes, PANELS):
    ax.scatter(df[xcol], df[ycol], s=4, alpha=0.1, color="tab:blue")
    ax.plot([AXIS_LO, AXIS_HI], [AXIS_LO, AXIS_HI], "--", color="gray", lw=1)
    ax.set_xlim(AXIS_LO, AXIS_HI)
    ax.set_ylim(AXIS_LO, AXIS_HI)
    ax.set_xlabel("default (68 tissues)")
    ax.set_ylabel("onebrain")
    ax.set_title(title, fontsize=10)
    annotate_diag_stats(ax, df[xcol], df[ycol])

    dist = (df[xcol] - df[ycol]).abs() / np.sqrt(2)
    ax.legend(handles=[plt.Line2D([0], [0], color="none",
                                   label=f"dist from y=x: mean={dist.mean():.4f}, median={dist.median():.4f}")],
              loc="lower right", handlelength=0, handletextpad=0, fontsize=8)

fig.suptitle(f"free_weights tau: default vs onebrain, n={len(df)} genes")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)
