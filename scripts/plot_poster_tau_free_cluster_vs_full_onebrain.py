"""Poster version of plot_tau_free_vs_full.py's row-1 panels only (gtex_protein_coding,
free_cluster_weights_max, onebrain leg): left, tau_unweighted (full 68 tissues) vs tau
on the brain-combined reference; right, tau_free_weighted (full) vs the same reference.
dpi=600, no suptitle, panel titles simplified to "tau unweighted"/"tau weighted"
(+9pt), pearson r line dropped from both panels' stats keys, "Brain" decapitalized
to "brain" in the category legend, all 3 per-panel keys (category legend,
diag-stats key, "median dist. from y=x:" key) at a matched font size, axis
labels sized up (+9pt), and tick labels sized up (+4pt) for poster legibility."""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from plot_utils import display_label
from gls_tau import floor_tpm

IN_GCT = sys.argv[1]
IN_FREE_TAU = sys.argv[2]
IN_REF_TAU = sys.argv[3]
OUT_PNG = sys.argv[4]

LABEL = "onebrain"
CATEGORIES = {"brain": "tab:red", "other": "tab:blue"}
AXIS_LO, AXIS_HI = 0.0, 1.0  # tau on log10(TPM+1) is lower; 0.2 clipped ~8% of genes
KEY_FONTSIZE = 9 + 4 + 4
LEGEND_FONTSIZE = KEY_FONTSIZE
AXIS_LABEL_FONTSIZE = 10 + 5 + 4
TITLE_FONTSIZE = 10 + 5 + 4
TICK_LABEL_FONTSIZE = 10 + 4

gct = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(gct.columns[2:])
mat = floor_tpm(gct[tissues].values.astype(float))
argmax_tissue = np.array(tissues)[mat.argmax(axis=1)]

category = np.full(len(argmax_tissue), "other", dtype=object)
category[np.char.find(argmax_tissue.astype(str), "Brain") >= 0] = "brain"
cat_df = pd.DataFrame({"gene": gct["Name"].values, "category": category})

d = pd.read_csv(IN_FREE_TAU)[["gene", "tau_unweighted", "tau_free_weighted"]]
r = pd.read_csv(IN_REF_TAU)[["gene", "tau_unweighted", "tau_weighted_diag_reg_inv_corr"]] \
    .rename(columns={"tau_unweighted": "tau_ref", "tau_weighted_diag_reg_inv_corr": "tau_weighted_alt"})

df = d.merge(r, on="gene").merge(cat_df, on="gene").dropna()

PANELS = [
    ("tau_unweighted", "tau_ref", "tau unweighted"),
    ("tau_free_weighted", "tau_ref", "tau weighted"),
]


def annotate_stats(ax, x, y):
    rho, _ = spearmanr(x, y)
    frac_below = (x > y).mean()
    frac_above = (x < y).mean()
    lines = []
    lines.append(f"spearman ρ={rho:.3f}")
    lines.append(f"below line={frac_below:.1%}")
    lines.append(f"above line={frac_above:.1%}")
    ax.text(0.03, 0.97, "\n".join(lines), transform=ax.transAxes, va="top", ha="left",
            fontsize=KEY_FONTSIZE, bbox=dict(facecolor="white", alpha=0.7, edgecolor="none"))


fig, axes = plt.subplots(1, 2, figsize=(19, 9.5))
for ax, (xcol, ycol, title) in zip(axes, PANELS):
    for cat, color in CATEGORIES.items():
        sub = df[df["category"] == cat]
        alpha = 0.1 if cat == "other" else 0.3
        ax.scatter(sub[xcol], sub[ycol], s=6, alpha=alpha, color=color, label=f"{cat} (n={len(sub)})")
    ax.plot([AXIS_LO, AXIS_HI], [AXIS_LO, AXIS_HI], "--", color="gray", lw=1)
    ax.set_xlim(AXIS_LO, AXIS_HI)
    ax.set_ylim(AXIS_LO, AXIS_HI)
    ax.set_xlabel("full (68 tissues)", fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel(display_label(LABEL), fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_title(title, fontsize=TITLE_FONTSIZE)
    ax.tick_params(labelsize=TICK_LABEL_FONTSIZE)
    annotate_stats(ax, df[xcol], df[ycol])

    cat_legend = ax.legend(loc="center left", fontsize=LEGEND_FONTSIZE, markerscale=3)
    ax.add_artist(cat_legend)

    dist = (df[xcol] - df[ycol]).abs() / np.sqrt(2)
    ax.legend(handles=[plt.Line2D([0], [0], color="none", label=f"median dist. from y=x: {dist.median():.4f}")],
              loc="lower right", handlelength=0, handletextpad=0, fontsize=KEY_FONTSIZE)

fig.tight_layout()
fig.savefig(OUT_PNG, dpi=600)
plt.close(fig)

print(f"-> {OUT_PNG}", file=sys.stderr)
