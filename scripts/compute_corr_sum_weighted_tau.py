"""Weighted tau using row-sum redundancy weights (see gls_tau.corr_sum_weights):
w_t = 1/(1 + sum_{t'!=t} r_tt'^POWER), from R = tissue x tissue Pearson correlation of
log10(TPM+1). POWER=1 ("corr_sum") uses raw correlation; POWER=2 ("corr_sq_sum") uses
r^2, so strong negative and positive correlations are penalized equally. Contrasted
with diag_reg_inv_corr (matrix inversion, needs diagonal regularization to avoid
negative weights) and margin (skips the correlation matrix entirely): this only needs
the row sums of R, no inversion, no regularization search. TPM<1 floored to 0 before
any computation. R is restricted to genes with TPM>=1 in at least one tissue (see
gls_tau.tissue_corr_expressed_only)."""
import sys
import numpy as np
import pandas as pd
from gls_tau import tissue_corr_expressed_only, corr_sum_weights, weighted_tau, floor_tpm

IN_GCT = sys.argv[1]
POWER = int(sys.argv[2])
OUT_TAU_CSV = sys.argv[3]
OUT_WEIGHTS_CSV = sys.argv[4]

METHOD_NAME = {1: "corr_sum", 2: "corr_sq_sum"}[POWER]

df = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(df.columns[2:])
genes = df["Name"].values
mat = df[tissues].values.astype(float)  # genes x tissues
mat = floor_tpm(mat)  # TPM<1 -> 0

R = tissue_corr_expressed_only(mat)
w = corr_sum_weights(R, power=POWER)

out = pd.DataFrame({
    "gene": genes,
    "tau_unweighted": weighted_tau(mat, np.ones(len(tissues))),
    f"tau_weighted_{METHOD_NAME}": weighted_tau(mat, w),
})
out.to_csv(OUT_TAU_CSV, index=False)

weights_out = pd.DataFrame({
    "tissue": tissues,
    f"weight_{METHOD_NAME}": w,
})
weights_out.to_csv(OUT_WEIGHTS_CSV, index=False)

print(f"genes={len(out)} tissues={len(tissues)} method={METHOD_NAME}", file=sys.stderr)
