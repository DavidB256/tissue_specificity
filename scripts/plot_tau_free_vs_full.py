"""2-panel comparison of tau on the full tissue set vs a reduced-resolution reference
(LABEL, e.g. "onebrain", "oneislets" or "hierclust0.1" -- unweighted tau on that
collapsed-tissue GCT, used as the optimization target, see optimize_free_weights.py):
(1) tau_unweighted vs tau_unweighted, (2) tau_free_weighted (full) vs tau_unweighted
(reference). Works for any source (gene or transcript rows). Points colored by gene category (by
max-expressing tissue, full 68-tissue GCT) vs "other", to show how the free_weights
correction differs across degrees of tissue redundancy (these tissue groups hold many
highly correlated subtissues). Category set depends on LABEL: onebrain uses Brain vs
other only; all other labels (e.g. hierclust0.1) use Brain/Kidney/Liver vs other."""
import sys
import gzip
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
IN_REF_GCT = sys.argv[6] if len(sys.argv) > 6 else None  # reduced-resolution GCT, for the tissue count

gct = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(gct.columns[2:])
if not any("Brain" in t for t in tissues):
    CATEGORIES = {"other": "tab:blue"}
elif LABEL == "onebrain":
    CATEGORIES = {"Brain": "tab:red", "other": "tab:blue"}
else:
    CATEGORIES = {"Brain": "tab:red", "Kidney": "tab:green", "Liver": "tab:orange", "other": "tab:blue"}
mat = floor_tpm(gct[tissues].values.astype(float))
argmax_tissue = np.array(tissues)[mat.argmax(axis=1)]

category = np.full(len(argmax_tissue), "other", dtype=object)
for cat in [c for c in CATEGORIES if c != "other"]:
    category[np.char.find(argmax_tissue.astype(str), cat) >= 0] = cat
cat_df = pd.DataFrame({"gene": gct["Name"].values, "category": category})

d = pd.read_csv(IN_FREE_TAU)[["gene", "tau_unweighted", "tau_free_weighted"]]
r = pd.read_csv(IN_REF_TAU)[["gene", "tau_unweighted"]].rename(columns={"tau_unweighted": "tau_ref"})

df = d.merge(r, on="gene").merge(cat_df, on="gene").dropna()

if IN_REF_GCT:
    UNIT = "cell types" if LABEL == "oneislets" else "tissues"
    MERGED_DESC = {"onebrain": "brain tissues merged", "oneislets": "islet cell types merged"}.get(LABEL, f"{LABEL} merged")
    n_ref = int(gzip.open(IN_REF_GCT, "rt").readlines()[1].split("\t")[1])
    REF_DESC = f"{MERGED_DESC}: {n_ref} {UNIT}"
    Y_LABEL = f"tau ({REF_DESC})"
    X_LABELS = [f"tau (full: {len(tissues)} {UNIT})", f"tau_weighted (full: {len(tissues)} {UNIT})"]
else:
    REF_DESC = display_label(LABEL)
    Y_LABEL = display_label(LABEL)
    X_LABELS = ["full", "full"]
PANELS = [
    ("tau_unweighted", "tau_ref", "tau_unweighted", X_LABELS[0], Y_LABEL),
    ("tau_free_weighted", "tau_ref", "tau_free_weighted", X_LABELS[1], Y_LABEL),
]

AXIS_LO, AXIS_HI = 0.0, 1.0  # tau on log10(TPM+1) is lower; 0.2 clipped ~8% of genes

fig, axes = plt.subplots(1, 2, figsize=(12, 6))
for ax, (xcol, ycol, title, xlabel, ylabel) in zip(axes, PANELS):
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
                                   label=f"median dist. from y=x: {dist.median():.4f}")],
              loc="lower right", handlelength=0, handletextpad=0, fontsize=8)

fig.suptitle(f"{OUT_PNG.split('/')[-4]} tau: full vs {REF_DESC}, n={len(df)} genes")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)
