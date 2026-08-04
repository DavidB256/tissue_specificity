"""PCA of tissues in one model's fitted (unconstrained) parameter space.
Plots PC1 vs PC2 beside PC1 vs PC3, every tissue labeled."""
import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mixture_core import theta_from_fit, canonicalize_by_weight, unflatten_fit
from model_registry import MODELS

MODEL = sys.argv[1]
IN_FITS = sys.argv[2]
OUT_PNG = sys.argv[3]
COMP_TYPES = MODELS[MODEL]

fits = pd.read_csv(IN_FITS)
tissues = fits["tissue"].tolist()
X = np.stack([
    theta_from_fit(COMP_TYPES, *canonicalize_by_weight(
        COMP_TYPES, *unflatten_fit(COMP_TYPES, row)))
    for _, row in fits.iterrows()
])

mean = X.mean(axis=0)
std = X.std(axis=0, ddof=1)
std[std < 1e-8] = 1.0
Xz = (X - mean) / std

U, S, Vt = np.linalg.svd(Xz, full_matrices=False)
scores = U * S
var_ratio = (S ** 2) / (S ** 2).sum()

fig, axes = plt.subplots(1, 2, figsize=(17, 8))
for ax, (pcx, pcy) in zip(axes, [(0, 1), (0, 2)]):
    ax.scatter(scores[:, pcx], scores[:, pcy], s=12, alpha=0.7, color="steelblue")
    for i, t in enumerate(tissues):
        ax.annotate(t, (scores[i, pcx], scores[i, pcy]), fontsize=4, alpha=0.85,
                    xytext=(2, 2), textcoords="offset points")
    ax.set_xlabel(f"PC{pcx + 1} ({var_ratio[pcx] * 100:.1f}% var)")
    ax.set_ylabel(f"PC{pcy + 1} ({var_ratio[pcy] * 100:.1f}% var)")
    ax.axhline(0, color="gray", lw=0.5)
    ax.axvline(0, color="gray", lw=0.5)
fig.suptitle(f"model={MODEL}: tissue PCA in fitted parameter space (k={X.shape[1]})", fontsize=13)
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)
