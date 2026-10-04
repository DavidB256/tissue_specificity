"""free_weights analog of plot_tau_cluster_vs_full.py: 4-panel comparison of tau on
full (68 tissues) vs a reduced-resolution reference (LABEL, e.g. "onebrain" or
"hierclust0.1" -- unweighted tau on that collapsed-tissue GCT, used as the
optimization target, see optimize_free_weights.py). Row 1: (1) tau_unweighted vs
tau_unweighted, (2) tau_free_weighted (full) vs tau_unweighted (reference). Row 2:
(3) tau vs tau_weighted on full tissues, (4) tau vs tau_weighted on the reference's
collapsed tissues (tau_weighted_diag_reg_inv_corr). Points colored by gene category (by
max-expressing tissue, full 68-tissue GCT) vs "other", to show how the free_weights
correction differs across degrees of tissue redundancy (these tissue groups hold many
highly correlated subtissues). Category set depends on LABEL: onebrain uses Brain vs
other only; all other labels (e.g. hierclust0.1) use Brain/Kidney/Liver vs other."""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plot_utils import annotate_diag_stats, display_label
from gls_tau import floor_tpm

IN_GCT = sys.argv[1]
IN_FREE_TAU = sys.argv[2]
IN_REF_TAU = sys.argv[3]
OUT_PNG = sys.argv[4]
LABEL = sys.argv[5] if len(sys.argv) > 5 else "reference"

if LABEL == "onebrain":
    CATEGORIES = {"Brain": "tab:red", "other": "tab:blue"}
else:
    CATEGORIES = {"Brain": "tab:red", "Kidney": "tab:green", "Liver": "tab:orange", "other": "tab:blue"}

gct = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(gct.columns[2:])
mat = floor_tpm(gct[tissues].values.astype(float))
argmax_tissue = np.array(tissues)[mat.argmax(axis=1)]

category = np.full(len(argmax_tissue), "other", dtype=object)
for cat in [c for c in CATEGORIES if c != "other"]:
    category[np.char.find(argmax_tissue.astype(str), cat) >= 0] = cat
cat_df = pd.DataFrame({"gene": gct["Name"].values, "category": category})

d = pd.read_csv(IN_FREE_TAU)[["gene", "tau_unweighted", "tau_free_weighted"]]
r = pd.read_csv(IN_REF_TAU)[["gene", "tau_unweighted", "tau_weighted_diag_reg_inv_corr"]] \
    .rename(columns={"tau_unweighted": "tau_ref", "tau_weighted_diag_reg_inv_corr": "tau_weighted_alt"})

df = d.merge(r, on="gene").merge(cat_df, on="gene").dropna()

ROW1_PANELS = [
    ("tau_unweighted", "tau_ref", "tau_unweighted", "full (68 tissues)", display_label(LABEL)),
    ("tau_free_weighted", "tau_ref", "tau_free_weighted", "full (68 tissues)", display_label(LABEL)),
]
ROW2_PANELS = [
    ("tau_unweighted", "tau_free_weighted", "full: tau vs tau_weighted",
     "tau (full)", "tau_weighted (full)"),
    ("tau_ref", "tau_weighted_alt", f"{display_label(LABEL)}: tau vs tau_weighted",
     f"tau ({display_label(LABEL)})", f"tau_weighted ({display_label(LABEL)})"),
]

AXIS_LO, AXIS_HI = 0.2, 1.0

fig, axes = plt.subplots(2, 2, figsize=(12, 12))
for row_axes, panels in zip(axes, [ROW1_PANELS, ROW2_PANELS]):
    for ax, (xcol, ycol, title, xlabel, ylabel) in zip(row_axes, panels):
        for cat, color in CATEGORIES.items():
            sub = df[df["category"] == cat]
            alpha = 0.1 if cat == "other" else 0.3
            ax.scatter(sub[xcol], sub[ycol], s=4, alpha=alpha, color=color, label=f"{cat} (n={len(sub)})")
        ax.plot([AXIS_LO, AXIS_HI], [AXIS_LO, AXIS_HI], "--", color="gray", lw=1)
        ax.set_xlim(AXIS_LO, AXIS_HI)
        ax.set_ylim(AXIS_LO, AXIS_HI)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title, fontsize=10)
        annotate_diag_stats(ax, df[xcol], df[ycol])

        cat_legend = ax.legend(loc="center left", fontsize=7, markerscale=3)
        ax.add_artist(cat_legend)

        dist = (df[xcol] - df[ycol]).abs() / np.sqrt(2)
        ax.legend(handles=[plt.Line2D([0], [0], color="none",
                                       label=f"dist from y=x: mean={dist.mean():.4f}, median={dist.median():.4f}")],
                  loc="lower right", handlelength=0, handletextpad=0, fontsize=8)

fig.suptitle(f"free_weights tau: full vs {display_label(LABEL)}, n={len(df)} genes")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)
