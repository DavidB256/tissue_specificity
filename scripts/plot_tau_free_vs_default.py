"""free_weights analog of plot_tau_cluster_vs_default.py: 2-panel comparison of tau on
default (68 tissues) vs a reduced-resolution reference (LABEL, e.g. "onebrain" or
"simplified_0.1" -- unweighted tau on that collapsed-tissue GCT, used as the
optimization target, see optimize_free_weights.py): (1) tau_unweighted vs tau_unweighted,
(2) tau_free_weighted (default) vs tau_unweighted (reference). Points colored by gene
category: Brain/Kidney/Liver (by max-expressing tissue, default 68-tissue GCT) vs
"other", to show how the free_weights correction differs across degrees of tissue
redundancy (Brain/Kidney/Liver hold many highly correlated subtissues)."""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plot_utils import annotate_diag_stats
from gls_tau import floor_tpm

IN_GCT = sys.argv[1]
IN_FREE_TAU = sys.argv[2]
IN_REF_TAU = sys.argv[3]
OUT_PNG = sys.argv[4]
LABEL = sys.argv[5] if len(sys.argv) > 5 else "reference"

CATEGORIES = {"Brain": "tab:red", "Kidney": "tab:green", "Liver": "tab:orange", "other": "tab:blue"}

gct = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(gct.columns[2:])
mat = floor_tpm(gct[tissues].values.astype(float))
argmax_tissue = np.array(tissues)[mat.argmax(axis=1)]

category = np.full(len(argmax_tissue), "other", dtype=object)
for cat in ["Brain", "Kidney", "Liver"]:
    category[np.char.find(argmax_tissue.astype(str), cat) >= 0] = cat
cat_df = pd.DataFrame({"gene": gct["Name"].values, "category": category})

d = pd.read_csv(IN_FREE_TAU)[["gene", "tau_unweighted", "tau_free_weighted"]]
r = pd.read_csv(IN_REF_TAU)[["gene", "tau_unweighted"]].rename(columns={"tau_unweighted": "tau_ref"})

df = d.merge(r, on="gene").merge(cat_df, on="gene").dropna()

PANELS = [
    ("tau_unweighted", "tau_ref", "tau_unweighted"),
    ("tau_free_weighted", "tau_ref", "tau_free_weighted"),
]

AXIS_LO, AXIS_HI = 0.2, 1.0

fig, axes = plt.subplots(1, 2, figsize=(12, 6))
for ax, (xcol, ycol, title) in zip(axes, PANELS):
    for cat, color in CATEGORIES.items():
        sub = df[df["category"] == cat]
        alpha = 0.1 if cat == "other" else 0.3
        ax.scatter(sub[xcol], sub[ycol], s=4, alpha=alpha, color=color, label=f"{cat} (n={len(sub)})")
    ax.plot([AXIS_LO, AXIS_HI], [AXIS_LO, AXIS_HI], "--", color="gray", lw=1)
    ax.set_xlim(AXIS_LO, AXIS_HI)
    ax.set_ylim(AXIS_LO, AXIS_HI)
    ax.set_xlabel("default (68 tissues)")
    ax.set_ylabel(LABEL)
    ax.set_title(title, fontsize=10)
    annotate_diag_stats(ax, df[xcol], df[ycol])

    cat_legend = ax.legend(loc="center left", fontsize=7, markerscale=3)
    ax.add_artist(cat_legend)

    dist = (df[xcol] - df[ycol]).abs() / np.sqrt(2)
    ax.legend(handles=[plt.Line2D([0], [0], color="none",
                                   label=f"dist from y=x: mean={dist.mean():.4f}, median={dist.median():.4f}")],
              loc="lower right", handlelength=0, handletextpad=0, fontsize=8)

fig.suptitle(f"free_weights tau: default vs {LABEL}, n={len(df)} genes")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)
