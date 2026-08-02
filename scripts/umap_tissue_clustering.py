"""UMAP of tissues in one model's fitted (unconstrained) parameter space.
Same standardized theta feature space as pca_tissue_clustering.py."""
import sys
import json
import pandas as pd
import numpy as np
import umap
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mixture_core import theta_from_fit, canonicalize_by_weight
from model_registry import MODELS

MODEL = sys.argv[1]
IN_FITS = sys.argv[2]
OUT_PNG = sys.argv[3]
COMP_TYPES = MODELS[MODEL]
N_NEIGHBORS = 15
MIN_DIST = 0.1
SEED = 42

fits = pd.read_csv(IN_FITS)
tissues = fits["tissue"].tolist()
X = np.stack([
    theta_from_fit(COMP_TYPES, *canonicalize_by_weight(
        COMP_TYPES, json.loads(row.weights), json.loads(row.params)))
    for _, row in fits.iterrows()
])

mean = X.mean(axis=0)
std = X.std(axis=0, ddof=1)
std[std < 1e-8] = 1.0
Xz = (X - mean) / std

reducer = umap.UMAP(n_components=2, n_neighbors=N_NEIGHBORS, min_dist=MIN_DIST, random_state=SEED)
embedding = reducer.fit_transform(Xz)

fig, ax = plt.subplots(figsize=(10, 9))
ax.scatter(embedding[:, 0], embedding[:, 1], s=14, alpha=0.7, color="steelblue")
for i, t in enumerate(tissues):
    ax.annotate(t, (embedding[i, 0], embedding[i, 1]), fontsize=4, alpha=0.85,
                xytext=(2, 2), textcoords="offset points")
ax.set_xlabel("UMAP1")
ax.set_ylabel("UMAP2")
ax.set_title(
    f"model={MODEL}: tissue UMAP in fitted parameter space "
    f"(k={X.shape[1]}, n={len(tissues)}, n_neighbors={N_NEIGHBORS}, min_dist={MIN_DIST})",
    fontsize=10,
)
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)
