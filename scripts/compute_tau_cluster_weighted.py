"""Weighted tau using cluster-size-based tissue weights instead of GLS: cut the
existing hierarchical clustering (tissue_correlation_linkage.csv, distance=1-corr,
average linkage -- see compute_tissue_correlation.py) at DIST_THRESHOLD (default 0.1,
i.e. tissues merge into the same cluster once corr>=0.9), then weight_i = 1/cluster_size
so every cluster's tissues sum to exactly 1 -- a cluster of k redundant tissues
contributes the same total weight as a single independent tissue, regardless of k. This
avoids GLS/ridge's global-lambda problem, where shrinkage needed to fix unrelated
negative weights elsewhere in R over-shrinks tight, unrelated clusters back toward
uniform (see cluster_weights design discussion)."""
import sys
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster
from gls_tau import weighted_tau, floor_tpm

IN_GCT = sys.argv[1]
IN_LINKAGE_CSV = sys.argv[2]
OUT_CLUSTERS_CSV = sys.argv[3]
OUT_TAU_CSV = sys.argv[4]
DIST_THRESHOLD = float(sys.argv[5]) if len(sys.argv) > 5 else 0.1

df = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(df.columns[2:])
genes = df["Name"].values
mat = floor_tpm(df[tissues].values.astype(float))

Z = pd.read_csv(IN_LINKAGE_CSV)[["cluster1", "cluster2", "distance", "n_obs"]].values
cluster_id = fcluster(Z, t=DIST_THRESHOLD, criterion="distance")

cluster_size = pd.Series(cluster_id).value_counts()
weights = np.array([1.0 / cluster_size[c] for c in cluster_id])

clusters_out = pd.DataFrame({
    "tissue": tissues,
    "cluster_id": cluster_id,
    "cluster_size": [cluster_size[c] for c in cluster_id],
    "weight": weights,
})
clusters_out.to_csv(OUT_CLUSTERS_CSV, index=False)

n_multi = int((clusters_out["cluster_size"] > 1).sum())
print(f"tissues={len(tissues)} clusters={cluster_id.max()} dist_threshold={DIST_THRESHOLD} "
      f"tissues_in_multi_clusters={n_multi}", file=sys.stderr)

out = pd.DataFrame({
    "gene": genes,
    "tau_unweighted": weighted_tau(mat, np.ones(len(tissues))),
    "tau_weighted_cluster": weighted_tau(mat, weights),
})
out.to_csv(OUT_TAU_CSV, index=False)

print(f"genes={len(out)} tissues={len(tissues)}", file=sys.stderr)
