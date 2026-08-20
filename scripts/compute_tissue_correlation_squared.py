"""Elementwise (Hadamard) square of the standard tissue x tissue correlation matrix
(gene_specificity/tissue_correlation.csv, from compute_tissue_correlation.py) --
R*R stays PSD (Schur product theorem) and keeps diagonal=1, up-weighting high |corr|
relative to moderate corr. Since all entries here are non-negative, squaring doesn't
flip sign. Clustering distance = 1 - R**2, average linkage, recomputed on the squared
matrix."""
import sys
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform

IN_CORR_CSV = sys.argv[1]
OUT_CORR_CSV = sys.argv[2]
OUT_LINKAGE_CSV = sys.argv[3]

corr_df = pd.read_csv(IN_CORR_CSV, index_col=0)
tissues = list(corr_df.columns)
R = corr_df.values
assert (R >= 0).all(), "squaring flips sign of negative correlations -- not valid here"

R2 = R ** 2
pd.DataFrame(R2, index=tissues, columns=tissues).to_csv(OUT_CORR_CSV)

dist = 1 - R2
np.fill_diagonal(dist, 0.0)
dist = (dist + dist.T) / 2
Z = linkage(squareform(dist, checks=False), method="average")
pd.DataFrame(Z, columns=["cluster1", "cluster2", "distance", "n_obs"]).to_csv(OUT_LINKAGE_CSV, index=False)

print(f"tissues={len(tissues)}", file=sys.stderr)
