"""PCA of tissues in raw expression space (log10(TPM+1), TPM<1->0, genes with TPM>=1
in at least one tissue -- same gene set as gls_tau.tissue_corr_expressed_only, so this
PCA lives in the same space the clustering (compute_tissue_correlation.py ->
compute_tau_cluster_weighted.py) was built from), colored by cluster_id from
tissue_clusters.csv. Centered-only (not standardized), same convention as
expr_vs_model_pca.py. Plots PC1 vs PC2 beside PC1 vs PC3, every tissue labeled."""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from gls_tau import floor_tpm

IN_GCT = sys.argv[1]
IN_CLUSTERS_CSV = sys.argv[2]
OUT_PNG = sys.argv[3]

gct = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(gct.columns[2:])
mat = floor_tpm(gct[tissues].values.astype(float))
expressed = mat.max(axis=1) > 0
log_mat = np.log10(mat[expressed] + 1.0)  # genes x tissues

X = log_mat.T  # tissues x genes
Xc = X - X.mean(axis=0)

U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
scores = U * S
var_ratio = (S ** 2) / (S ** 2).sum()

clusters = pd.read_csv(IN_CLUSTERS_CSV).set_index("tissue")
cluster_id = clusters.loc[tissues, "cluster_id"].values
cluster_size = clusters.loc[tissues, "cluster_size"].values
uniq_clusters = sorted(set(cluster_id))
cmap = plt.get_cmap("tab20")
color_of = {c: cmap(i % 20) for i, c in enumerate(uniq_clusters)}
colors = [color_of[c] for c in cluster_id]

fig, axes = plt.subplots(1, 2, figsize=(19, 8))
for ax, (pcx, pcy) in zip(axes, [(0, 1), (0, 2)]):
    for c in uniq_clusters:
        mask = cluster_id == c
        n = cluster_size[mask][0]
        ax.scatter(scores[mask, pcx], scores[mask, pcy], s=20, alpha=0.85,
                   color=color_of[c], label=f"cluster {c} (n={n})")
    for i, t in enumerate(tissues):
        ax.annotate(t, (scores[i, pcx], scores[i, pcy]), fontsize=4, alpha=0.85,
                    xytext=(2, 2), textcoords="offset points")
    ax.set_xlabel(f"PC{pcx + 1} ({var_ratio[pcx] * 100:.1f}% var)")
    ax.set_ylabel(f"PC{pcy + 1} ({var_ratio[pcy] * 100:.1f}% var)")
    ax.axhline(0, color="gray", lw=0.5)
    ax.axvline(0, color="gray", lw=0.5)

axes[1].legend(fontsize=5, ncol=2, loc="upper left", bbox_to_anchor=(1.01, 1.0))
fig.suptitle(f"tissue PCA in expression space (log10 TPM, TPM<1->0), colored by cluster (dist<=0.1)", fontsize=13)
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)
