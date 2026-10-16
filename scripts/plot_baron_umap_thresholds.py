"""3-panel figure: UMAP of Baron et al. human pancreas single cells (protein-
coding TPM, compute_baron_cell_tpm.py), one panel per TPM-floor choice --
never floor, TPM<1->0, TPM<10->0 -- log10(TPM+1) -> z-score -> PCA(30) -> UMAP,
same pipeline/seed across panels except for the floor step, isolating its
effect. Cell-type clusters are CALLED by KMeans (k=14, matching the number of
Baron-labeled cell types) on the PCA space, not read from Baron's own labels;
each panel's title reports the adjusted Rand index between the called KMeans
clusters and Baron's assigned_cluster labels, as a check of how well each
floor choice recovers the known cell types."""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
from umap import UMAP
from gls_tau import floor_tpm

IN_CELL_TPM_CSV = sys.argv[1]
OUT_PNG = sys.argv[2]

SEED = 0
N_PCS = 30
THRESHOLDS = [("never round down", 0.0), ("TPM<1 -> 0", 1.0), ("TPM<10 -> 0", 10.0)]

df = pd.read_csv(IN_CELL_TPM_CSV)
assigned_cluster = df["assigned_cluster"].to_numpy()
n_true_types = len(np.unique(assigned_cluster))
gene_cols = [c for c in df.columns if c not in ("assigned_cluster", "sample")]
raw_mat = df[gene_cols].to_numpy(dtype=float)
print(f"cells={len(df)} genes={len(gene_cols)} true_cell_types={n_true_types}", file=sys.stderr)

fig, axes = plt.subplots(1, 3, figsize=(19, 6.5))

for ax, (label, floor) in zip(axes, THRESHOLDS):
    mat = floor_tpm(raw_mat, floor=floor) if floor > 0 else raw_mat
    log_mat = np.log10(mat + 1.0)
    scaled = StandardScaler().fit_transform(log_mat)
    pcs = PCA(n_components=N_PCS, random_state=SEED).fit_transform(scaled)
    emb = UMAP(random_state=SEED).fit_transform(pcs)
    called = KMeans(n_clusters=n_true_types, random_state=SEED, n_init=10).fit_predict(pcs)
    ari = adjusted_rand_score(assigned_cluster, called)

    sca = ax.scatter(emb[:, 0], emb[:, 1], c=called, cmap="tab20", s=4, alpha=0.8)
    ax.set_title(f"{label}\nKMeans (k={n_true_types}) vs. Baron labels: ARI={ari:.3f}")
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    print(f"{label}: ari={ari:.4f}", file=sys.stderr)

fig.suptitle("Baron pancreas scRNA-seq: UMAP + called clusters by TPM floor choice")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)

print(f"-> {OUT_PNG}", file=sys.stderr)
