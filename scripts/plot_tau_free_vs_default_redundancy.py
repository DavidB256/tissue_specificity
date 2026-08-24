"""plot_tau_free_vs_default.py colored by gene category: Brain/Kidney/Liver (by
max-expressing tissue, default 68-tissue GCT) vs "other". Unlike the older _bkl figures,
"other" genes are kept (not dropped) so the coloring shows how the free_weights
correction differs across degrees of tissue redundancy. Panels compare default (68
tissues) against three one-tissue-collapsed scenarios -- onebrain, onekidney, oneliver
(see build_gct_onebrain.py / build_gct_onekidney.py / build_gct_oneliver.py) -- reusing
the single free_weights fit (optimized against onebrain only; no weights are recomputed
for onekidney/oneliver here)."""
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
IN_ONEBRAIN_TAU = sys.argv[3]
IN_ONEKIDNEY_TAU = sys.argv[4]
IN_ONELIVER_TAU = sys.argv[5]
OUT_PNG = sys.argv[6]

CATEGORIES = {"Brain": "tab:red", "Kidney": "tab:green", "Liver": "tab:orange", "other": "tab:blue"}
SCENARIOS = ["onebrain", "onekidney", "oneliver"]

gct = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(gct.columns[2:])
mat = floor_tpm(gct[tissues].values.astype(float))
argmax_tissue = np.array(tissues)[mat.argmax(axis=1)]

category = np.full(len(argmax_tissue), "other", dtype=object)
for cat in ["Brain", "Kidney", "Liver"]:
    category[np.char.find(argmax_tissue.astype(str), cat) >= 0] = cat
cat_df = pd.DataFrame({"gene": gct["Name"].values, "category": category})

d = pd.read_csv(IN_FREE_TAU)[["gene", "tau_unweighted", "tau_free_weighted"]]
refs = {
    "onebrain": pd.read_csv(IN_ONEBRAIN_TAU)[["gene", "tau_unweighted"]].rename(columns={"tau_unweighted": "tau_onebrain"}),
    "onekidney": pd.read_csv(IN_ONEKIDNEY_TAU)[["gene", "tau_unweighted"]].rename(columns={"tau_unweighted": "tau_onekidney"}),
    "oneliver": pd.read_csv(IN_ONELIVER_TAU)[["gene", "tau_unweighted"]].rename(columns={"tau_unweighted": "tau_oneliver"}),
}

df = d.merge(cat_df, on="gene")
for ref in refs.values():
    df = df.merge(ref, on="gene")
df = df.dropna()

ROWS = [("tau_unweighted", "tau_unweighted"), ("tau_free_weighted", "tau_free_weighted")]

AXIS_LO, AXIS_HI = 0.2, 1.0

fig, axes = plt.subplots(2, 3, figsize=(21, 14))
for row, (xcol, row_label) in enumerate(ROWS):
    for col, scenario in enumerate(SCENARIOS):
        ax = axes[row, col]
        ycol = f"tau_{scenario}"
        for cat, color in CATEGORIES.items():
            sub = df[df["category"] == cat]
            alpha = 0.1 if cat == "other" else 0.3
            ax.scatter(sub[xcol], sub[ycol], s=4, alpha=alpha, color=color, label=f"{cat} (n={len(sub)})")
        ax.plot([AXIS_LO, AXIS_HI], [AXIS_LO, AXIS_HI], "--", color="gray", lw=1)
        ax.set_xlim(AXIS_LO, AXIS_HI)
        ax.set_ylim(AXIS_LO, AXIS_HI)
        ax.set_xlabel("default (68 tissues)")
        ax.set_ylabel(scenario)
        ax.set_title(f"{row_label} vs {scenario}", fontsize=10)
        annotate_diag_stats(ax, df[xcol], df[ycol])

        if row == 0 and col == 0:
            cat_legend = ax.legend(loc="upper left", fontsize=7, markerscale=3)
            ax.add_artist(cat_legend)

        dist = (df[xcol] - df[ycol]).abs() / np.sqrt(2)
        ax.legend(handles=[plt.Line2D([0], [0], color="none",
                                       label=f"dist from y=x: mean={dist.mean():.4f}, median={dist.median():.4f}")],
                  loc="lower right", handlelength=0, handletextpad=0, fontsize=8)

fig.suptitle(f"free_weights tau: default vs onebrain/onekidney/oneliver, colored by tissue redundancy, n={len(df)} genes")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)
