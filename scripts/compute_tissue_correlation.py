"""Tissue x tissue correlation matrix (human GTEx only) + hierarchical clustering.
TPM<1 rounded down to 0 before computing correlation (matches TPM>=1 expression floor
used elsewhere in the pipeline). Correlation on log10(TPM+1) (see gls_tau.tissue_corr).
Clustering distance = 1 - correlation, average linkage."""
import sys
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform
from gls_tau import tissue_corr, floor_tpm

IN_GCT = sys.argv[1]
OUT_CORR_CSV = sys.argv[2]
OUT_LINKAGE_CSV = sys.argv[3]

df = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(df.columns[2:])
mat = floor_tpm(df[tissues].values.astype(float))

R = tissue_corr(mat)
pd.DataFrame(R, index=tissues, columns=tissues).to_csv(OUT_CORR_CSV)

dist = 1 - R
np.fill_diagonal(dist, 0.0)
dist = (dist + dist.T) / 2  # enforce exact symmetry against float error
Z = linkage(squareform(dist, checks=False), method="average")
pd.DataFrame(Z, columns=["cluster1", "cluster2", "distance", "n_obs"]).to_csv(OUT_LINKAGE_CSV, index=False)

print(f"tissues={len(tissues)}", file=sys.stderr)
