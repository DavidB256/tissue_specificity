"""Cluster-size downweighting: alternative to GLS/ridge tissue weighting. Each tissue's
weight is 1/block_size, from the tissue_clusters_{k}.csv correlation-dendrogram cut (see
identify_tissue_clusters.py) -- tissues in bigger (more redundant) blocks get
downweighted proportionally to block size. No matrix inversion or regularization
involved, so weights are always positive. weighted_tau's w.sum()-w[argmax]
normalization makes the absolute scale of w irrelevant, so 1/block_size needs no
further normalization."""
import sys
import numpy as np
import pandas as pd
from gls_tau import weighted_tau, floor_tpm

IN_GCT = sys.argv[1]
IN_CLUSTERS = sys.argv[2]
OUT_TAU_CSV = sys.argv[3]
OUT_WEIGHTS_CSV = sys.argv[4]

df = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(df.columns[2:])
genes = df["Name"].values
mat = floor_tpm(df[tissues].values.astype(float))

clusters = pd.read_csv(IN_CLUSTERS).set_index("tissue").loc[tissues]
block_size = clusters["block"].map(clusters["block"].value_counts())
w = (1.0 / block_size).values

out = pd.DataFrame({
    "gene": genes,
    "tau_unweighted": weighted_tau(mat, np.ones(len(tissues))),
    "tau_clustersize": weighted_tau(mat, w),
})
out.to_csv(OUT_TAU_CSV, index=False)

weights_out = pd.DataFrame({
    "tissue": tissues,
    "block": clusters["block"].values,
    "block_size": block_size.values,
    "weight_clustersize": w,
})
weights_out.to_csv(OUT_WEIGHTS_CSV, index=False)

print(f"genes={len(out)} tissues={len(tissues)} n_blocks={clusters['block'].nunique()}", file=sys.stderr)
