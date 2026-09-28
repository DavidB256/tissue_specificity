"""2D histograms (20x20 bins) of tau vs frac_brain (fraction of a gene's total
floor_tpm'd linear TPM sitting in Brain_* tissues, see
list_neural_threshold_crossing_genes.py), one panel for tau_unweighted and one for
tau_free_weighted, same gene population (all genes with tau defined) side by side.
Args: IN_GCT IN_TAU_CSV OUT_PNG"""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from gls_tau import floor_tpm

IN_GCT = sys.argv[1]
IN_TAU_CSV = sys.argv[2]
OUT_PNG = sys.argv[3]

gct = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(gct.columns[2:])
mat = floor_tpm(gct[tissues].values.astype(float))
brain_idx = np.array([i for i, t in enumerate(tissues) if "Brain" in t])
total = mat.sum(axis=1)
brain_total = mat[:, brain_idx].sum(axis=1)
with np.errstate(invalid="ignore", divide="ignore"):
    frac_brain = np.where(total > 0, brain_total / total, np.nan)

genes = pd.DataFrame({"gene": gct["Name"].values, "frac_brain": frac_brain})
tau = pd.read_csv(IN_TAU_CSV)[["gene", "tau_unweighted", "tau_free_weighted"]]
df = tau.merge(genes, on="gene").dropna()

fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
panels = [
    ("tau_unweighted", "tau_unweighted vs frac_brain"),
    ("tau_free_weighted", "tau_weighted (free_cluster) vs frac_brain"),
]
for ax, (col, title) in zip(axes, panels):
    h = ax.hist2d(df["frac_brain"], df[col], bins=20, range=[[0, 1], [0, 1]], cmap="viridis")
    fig.colorbar(h[3], ax=ax, label="gene count")
    ax.set_xlabel("frac_brain (linear TPM, floor_tpm'd)")
    ax.set_ylabel(col)
    ax.set_title(f"{title} (n={len(df)} genes)")

fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)

print(f"genes={len(df)} -> {OUT_PNG}", file=sys.stderr)
