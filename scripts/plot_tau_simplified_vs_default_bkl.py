"""plot_tau_simplified_vs_default.py restricted to genes whose max-expressing tissue
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
IN_DEFAULT_TAU_SQUARED = sys.argv[3]
IN_SIMPLIFIED_TAU = sys.argv[4]
IN_SIMPLIFIED_TAU_SQUARED = sys.argv[5]
OUT_PNG = sys.argv[6]

CATEGORIES = {"Brain": "tab:red", "Kidney": "tab:green", "Liver": "tab:orange"}

gct = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(gct.columns[2:])
mat = floor_tpm(gct[tissues].values.astype(float))
argmax_tissue = np.array(tissues)[mat.argmax(axis=1)]

category = np.full(len(argmax_tissue), None, dtype=object)
for cat in CATEGORIES:
    category[np.char.find(argmax_tissue.astype(str), cat) >= 0] = cat
cat_df = pd.DataFrame({"gene": gct["Name"].values, "category": category})

d = pd.read_csv(IN_DEFAULT_TAU)[["gene", "tau_unweighted", "tau_weighted_gls_ridge"]]
d_sq = pd.read_csv(IN_DEFAULT_TAU_SQUARED)[["gene", "tau_weighted_gls_ridge"]].rename(
    columns={"tau_weighted_gls_ridge": "tau_weighted_gls_ridge_sq"})
s = pd.read_csv(IN_SIMPLIFIED_TAU)[["gene", "tau_unweighted", "tau_weighted_gls_ridge"]].rename(
    columns={"tau_unweighted": "tau_unweighted_simp", "tau_weighted_gls_ridge": "tau_weighted_gls_ridge_simp"})
s_sq = pd.read_csv(IN_SIMPLIFIED_TAU_SQUARED)[["gene", "tau_weighted_gls_ridge"]].rename(
    columns={"tau_weighted_gls_ridge": "tau_weighted_gls_ridge_simp_sq"})

df = d.merge(d_sq, on="gene").merge(s, on="gene").merge(s_sq, on="gene").merge(cat_df, on="gene").dropna()
df = df[df["category"].notna()]

PANELS = [
    ("tau_unweighted", "tau_unweighted_simp", "tau_unweighted: default vs simplified tissues"),
    ("tau_weighted_gls_ridge", "tau_weighted_gls_ridge_simp",
     "tau_weighted_gls_ridge (default R): default vs simplified tissues"),
    ("tau_weighted_gls_ridge_sq", "tau_weighted_gls_ridge_simp_sq",
     "tau_weighted_gls_ridge (squared R): default vs simplified tissues"),
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
    ax.set_ylabel(f"{ycol} (simplified tissues)")
    ax.set_title(title, fontsize=10)
    annotate_diag_stats(ax, df[xcol], df[ycol])

    cat_legend = ax.legend(loc="upper left", fontsize=7, markerscale=3)
    ax.add_artist(cat_legend)

    dist = (df[xcol] - df[ycol]).abs() / np.sqrt(2)
    ax.legend(handles=[plt.Line2D([0], [0], color="none",
                                   label=f"dist from y=x: mean={dist.mean():.4f}, median={dist.median():.4f}")],
              loc="lower right", handlelength=0, handletextpad=0, fontsize=8)

fig.suptitle(f"tau: default vs simplified tissues, max-expr in Brain/Kidney/Liver, n={len(df)} genes")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)
