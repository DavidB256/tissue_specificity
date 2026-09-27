"""Weighted tau using margin weights (see gls_tau.margin_weights): each tissue's
weight is the sum, over genes where it's the top expressor, of the margin to the
runner-up tissue. Redundant tissues (correlated groups) rarely win by a wide margin
and so accumulate little weight; distinctive tissues win decisively and accumulate
more. Contrasted with free/free_cluster (DE-optimized against a reference tau) and
diag_reg_inv_corr (closed-form from the tissue x tissue correlation matrix): margin
weights are closed-form directly from the expression matrix itself, no correlation
matrix, no reference tau, no regularization -- and never negative by construction.
TPM<1 floored to 0 before any computation."""
import sys
import numpy as np
import pandas as pd
from gls_tau import margin_weights, weighted_tau, floor_tpm

IN_GCT = sys.argv[1]
OUT_TAU_CSV = sys.argv[2]
OUT_WEIGHTS_CSV = sys.argv[3]

df = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(df.columns[2:])
genes = df["Name"].values
mat = df[tissues].values.astype(float)  # genes x tissues
mat = floor_tpm(mat)  # TPM<1 -> 0

w = margin_weights(mat)

out = pd.DataFrame({
    "gene": genes,
    "tau_unweighted": weighted_tau(mat, np.ones(len(tissues))),
    "tau_weighted_margin": weighted_tau(mat, w),
})
out.to_csv(OUT_TAU_CSV, index=False)

weights_out = pd.DataFrame({
    "tissue": tissues,
    "weight_margin": w,
})
weights_out.to_csv(OUT_WEIGHTS_CSV, index=False)

print(f"genes={len(out)} tissues={len(tissues)}", file=sys.stderr)
