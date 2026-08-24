"""Weighted tau using cluster-based tissue weights instead of GLS: cut the existing
hierarchical clustering (tissue_correlation_linkage.csv, distance=1-corr, average
linkage -- see compute_tissue_correlation.py) at DIST_THRESHOLD (default 0.1, i.e.
tissues merge into the same cluster once corr>=0.9). Two weighting schemes computed
side by side:
  weight        = 1/cluster_size -- every cluster's tissues sum to exactly 1, so a
                  cluster of k redundant tissues contributes the same total weight as
                  a single independent tissue regardless of k.
  weight_diameter = (1-d)/cluster_size + WEIGHT_CAP*d, where d = min(diameter/DIST_THRESHOLD, 1)
                  and diameter = max pairwise distance (1-corr) between any two tissues
                  in the cluster (from the actual correlation matrix, not just the
                  merge height -- average linkage can chain tissues into one cluster
                  even if some pairs inside it are less tightly correlated than the
                  cutoff). d=0 (tissues truly identical) -> full 1/k collapse; d=1
                  (diameter at or past the cutoff) -> weight WEIGHT_CAP (<1), never
                  fully independent -- clustering is never fully undone, even for the
                  loosest admitted clusters. Lets loose, chained clusters keep more of
                  their independent weight than tight ones, without ever discarding the
                  clustering correction entirely.
Both avoid GLS/ridge's global-lambda problem, where shrinkage needed to fix unrelated
negative weights elsewhere in R over-shrinks tight, unrelated clusters back toward
uniform (see cluster_weights design discussion)."""
import sys
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster
from gls_tau import weighted_tau, floor_tpm

IN_GCT = sys.argv[1]
IN_LINKAGE_CSV = sys.argv[2]
IN_CORR_CSV = sys.argv[3]
OUT_CLUSTERS_CSV = sys.argv[4]
OUT_TAU_CSV = sys.argv[5]
DIST_THRESHOLD = float(sys.argv[6]) if len(sys.argv) > 6 else 0.1
WEIGHT_CAP = float(sys.argv[7]) if len(sys.argv) > 7 else 0.5

df = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(df.columns[2:])
genes = df["Name"].values
mat = floor_tpm(df[tissues].values.astype(float))

Z = pd.read_csv(IN_LINKAGE_CSV)[["cluster1", "cluster2", "distance", "n_obs"]].values
cluster_id = fcluster(Z, t=DIST_THRESHOLD, criterion="distance")

R = pd.read_csv(IN_CORR_CSV, index_col=0).loc[tissues, tissues].values
dist_mat = 1 - R

cluster_size = pd.Series(cluster_id).value_counts()
weights = np.array([1.0 / cluster_size[c] for c in cluster_id])

diameter = np.zeros(len(tissues))
for c in cluster_size.index:
    idx = np.where(cluster_id == c)[0]
    if len(idx) > 1:
        diameter[idx] = dist_mat[np.ix_(idx, idx)].max()
d = np.minimum(diameter / DIST_THRESHOLD, 1.0)
weights_diameter = (1 - d) / np.array([cluster_size[c] for c in cluster_id]) + WEIGHT_CAP * d

clusters_out = pd.DataFrame({
    "tissue": tissues,
    "cluster_id": cluster_id,
    "cluster_size": [cluster_size[c] for c in cluster_id],
    "weight": weights,
    "diameter": diameter,
    "d": d,
    "weight_diameter": weights_diameter,
})
clusters_out.to_csv(OUT_CLUSTERS_CSV, index=False)

n_multi = int((clusters_out["cluster_size"] > 1).sum())
print(f"tissues={len(tissues)} clusters={cluster_id.max()} dist_threshold={DIST_THRESHOLD} "
      f"weight_cap={WEIGHT_CAP} tissues_in_multi_clusters={n_multi}", file=sys.stderr)

out = pd.DataFrame({
    "gene": genes,
    "tau_unweighted": weighted_tau(mat, np.ones(len(tissues))),
    "tau_weighted_cluster": weighted_tau(mat, weights),
    "tau_weighted_cluster_diameter": weighted_tau(mat, weights_diameter),
})
out.to_csv(OUT_TAU_CSV, index=False)

print(f"genes={len(out)} tissues={len(tissues)}", file=sys.stderr)
