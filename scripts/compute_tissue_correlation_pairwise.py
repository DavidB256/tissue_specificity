"""Pairwise-masked variant of compute_tissue_correlation.py: each tissue-pair's
correlation is restricted to genes with TPM>=1 (post floor) in at least one of the two
tissues, rather than one fixed gene set for the whole matrix. This violates the
single-consistent-sample-matrix assumption behind PSD, so R here is NOT guaranteed
PSD -- diagnostic use only, do not feed into gls_weights/ridge_lambda_search."""
import sys
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform
from gls_tau import floor_tpm

IN_GCT = sys.argv[1]
OUT_CORR_CSV = sys.argv[2]
OUT_LINKAGE_CSV = sys.argv[3]

df = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(df.columns[2:])
mat = floor_tpm(df[tissues].values.astype(float))
log_mat = np.log10(mat + 1.0)

n = len(tissues)
R = np.eye(n)
for i in range(n):
    for j in range(i + 1, n):
        mask = (mat[:, i] > 0) | (mat[:, j] > 0)
        r = np.corrcoef(log_mat[mask, i], log_mat[mask, j])[0, 1]
        R[i, j] = R[j, i] = r
    print(f"tissue {i + 1}/{n} done", file=sys.stderr)

pd.DataFrame(R, index=tissues, columns=tissues).to_csv(OUT_CORR_CSV)

dist = 1 - R
np.fill_diagonal(dist, 0.0)
dist = (dist + dist.T) / 2  # enforce exact symmetry against float error
Z = linkage(squareform(dist, checks=False), method="average")
pd.DataFrame(Z, columns=["cluster1", "cluster2", "distance", "n_obs"]).to_csv(OUT_LINKAGE_CSV, index=False)

print(f"tissues={n}", file=sys.stderr)
