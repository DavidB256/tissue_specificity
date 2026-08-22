"""plot_tau_cluster_vs_default.py restricted to genes whose max-expressing tissue
(default, unsimplified, 68 tissues) is Brain_*, Kidney_*, or Liver* -- colored by which
of the three."""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plot_utils import annotate_diag_stats
from gls_tau import floor_tpm

IN_GCT = sys.argv[1]
IN_DEFAULT_TAU = sys.argv[2]
IN_COLLAPSED_TAU = sys.argv[3]
OUT_PNG = sys.argv[4]
LABEL = sys.argv[5] if len(sys.argv) > 5 else "onebrain"

CATEGORIES = {"Brain": "tab:red", "Kidney": "tab:green", "Liver": "tab:orange"}

gct = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(gct.columns[2:])
mat = floor_tpm(gct[tissues].values.astype(float))
argmax_tissue = np.array(tissues)[mat.argmax(axis=1)]

category = np.full(len(argmax_tissue), None, dtype=object)
for cat in CATEGORIES:
    category[np.char.find(argmax_tissue.astype(str), cat) >= 0] = cat
cat_df = pd.DataFrame({"gene": gct["Name"].values, "category": category})

COLS = ["tau_unweighted", "tau_weighted_cluster", "tau_weighted_cluster_diameter"]
d = pd.read_csv(IN_DEFAULT_TAU)[["gene"] + COLS]
c = pd.read_csv(IN_COLLAPSED_TAU)[["gene"] + COLS].rename(columns={col: f"{col}_simp" for col in COLS})

df = d.merge(c, on="gene").merge(cat_df, on="gene").dropna()
df = df[df["category"].notna()]

PANELS = [
    ("tau_unweighted", "tau_unweighted_simp", f"tau_unweighted: default vs {LABEL} tissues"),
    ("tau_weighted_cluster", "tau_weighted_cluster_simp",
     f"tau_weighted_cluster (1/cluster_size): default vs {LABEL} tissues"),
    ("tau_weighted_cluster_diameter", "tau_weighted_cluster_diameter_simp",
     f"tau_weighted_cluster_diameter ((1-d)/k+d): default vs {LABEL} tissues"),
]

AXIS_LO, AXIS_HI = 0.2, 1.0

fig, axes = plt.subplots(1, 3, figsize=(18, 6))
for ax, (xcol, ycol, title) in zip(axes, PANELS):
    for cat, color in CATEGORIES.items():
        sub = df[df["category"] == cat]
        ax.scatter(sub[xcol], sub[ycol], s=4, alpha=0.1, color=color, label=f"{cat} (n={len(sub)})")
    ax.plot([AXIS_LO, AXIS_HI], [AXIS_LO, AXIS_HI], "--", color="gray", lw=1)
    ax.set_xlim(AXIS_LO, AXIS_HI)
    ax.set_ylim(AXIS_LO, AXIS_HI)
    ax.set_xlabel(f"{xcol} (default, 68 tissues)")
    ax.set_ylabel(f"{ycol} ({LABEL} tissues)")
    ax.set_title(title, fontsize=10)
    annotate_diag_stats(ax, df[xcol], df[ycol])

    cat_legend = ax.legend(loc="upper left", fontsize=7, markerscale=3)
    ax.add_artist(cat_legend)

    dist = (df[xcol] - df[ycol]).abs() / np.sqrt(2)
    ax.legend(handles=[plt.Line2D([0], [0], color="none",
                                   label=f"dist from y=x: mean={dist.mean():.4f}, median={dist.median():.4f}")],
              loc="lower right", handlelength=0, handletextpad=0, fontsize=8)

fig.suptitle(f"cluster-weighted tau: default vs {LABEL} tissues, max-expr in Brain/Kidney/Liver, n={len(df)} genes")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)
