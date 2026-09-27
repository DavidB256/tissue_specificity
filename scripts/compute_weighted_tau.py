"""Weighted tau: extends tau (Yanai et al. 2005) with diag_reg_inv_corr ("diagonal-
regularized inverse correlation") weights that account for cross-tissue correlation
(e.g. redundant brain subregions), so correlated tissues don't each pull full
independent weight in the specificity average. R = tissue x tissue Pearson
correlation of log10(TPM+1), computed across genes. Weights w solve R_lambda w = 1
(R_lambda = (1-lambda)*R + lambda*I) at the smallest lambda that makes every weight
positive (see gls_tau.diag_reg_lambda_search); w=1 for every tissue recovers
unweighted tau exactly when R=I, so no extra rescaling is needed. The unregularized
(lambda=0) solve is not reported on its own -- it routinely produces negative weights
under multicollinearity, which is not useful for this application. Ties for a gene's
max-expressing tissue go to the first occurrence (np.argmax default). TPM<1 floored to
0 before any computation. R is restricted to genes with TPM>=1 in at least one tissue
(see gls_tau.tissue_corr_expressed_only)."""
import sys
import numpy as np
import pandas as pd
from gls_tau import tissue_corr_expressed_only, diag_reg_lambda_search, weighted_tau, floor_tpm

IN_GCT = sys.argv[1]
OUT_TAU_CSV = sys.argv[2]
OUT_WEIGHTS_CSV = sys.argv[3]

df = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(df.columns[2:])
genes = df["Name"].values
mat = df[tissues].values.astype(float)  # genes x tissues
mat = floor_tpm(mat)  # TPM<1 -> 0

R = tissue_corr_expressed_only(mat)
lam_star, w_reg = diag_reg_lambda_search(R)
print(f"diag_reg lambda*={lam_star:.6g}", file=sys.stderr)

out = pd.DataFrame({
    "gene": genes,
    "tau_unweighted": weighted_tau(mat, np.ones(len(tissues))),
    "tau_weighted_diag_reg_inv_corr": weighted_tau(mat, w_reg),
})
out.to_csv(OUT_TAU_CSV, index=False)

weights_out = pd.DataFrame({
    "tissue": tissues,
    "weight_diag_reg_inv_corr": w_reg,
    "diag_reg_lambda": lam_star,
})
weights_out.to_csv(OUT_WEIGHTS_CSV, index=False)

print(f"genes={len(out)} tissues={len(tissues)} diag_reg_lambda={lam_star:.6g}", file=sys.stderr)
