"""Weighted tau: extends tau (Yanai et al. 2005) with GLS weights that account for
cross-tissue correlation (e.g. redundant brain subregions), so correlated tissues
don't each pull full independent weight in the specificity average. R = tissue x
tissue Pearson correlation of log10(TPM+1), computed across genes. GLS weights w
solve R w = 1 (proportional to R^-1 1); w=1 for every tissue recovers unweighted tau
exactly when R=I, so no extra rescaling is needed. Negative entries in w are possible
under multicollinearity and are reported, not clipped. diag_reg_inv_corr ("diagonal-
regularized inverse correlation") weights use R_lambda = (1-lambda)*R + lambda*I with
the smallest lambda that makes every weight positive (see gls_tau.diag_reg_lambda_search).
Ties for a gene's max-expressing tissue
go to the first occurrence (np.argmax default). TPM<1 floored to 0 before any
computation. R is restricted to genes with TPM>=1 in at least one tissue (see
gls_tau.tissue_corr_expressed_only). tau_weighted_gls (unregularized, unbounded below)
is left-censored at -1 to drop runaway outliers we don't care about."""
import sys
import numpy as np
import pandas as pd
from gls_tau import tissue_corr_expressed_only, gls_weights, diag_reg_lambda_search, weighted_tau, floor_tpm, censor_floor

IN_GCT = sys.argv[1]
OUT_TAU_CSV = sys.argv[2]
OUT_WEIGHTS_CSV = sys.argv[3]

df = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(df.columns[2:])
genes = df["Name"].values
mat = df[tissues].values.astype(float)  # genes x tissues
mat = floor_tpm(mat)  # TPM<1 -> 0

R = tissue_corr_expressed_only(mat)
w_unreg = gls_weights(R, 0.0)
neg_mask = w_unreg < 0
n_negative = int(neg_mask.sum())
if n_negative:
    print(f"negative unregularized GLS weights ({n_negative}/{len(tissues)}): "
          f"{[t for t, m in zip(tissues, neg_mask) if m]}", file=sys.stderr)

lam_star, w_reg = diag_reg_lambda_search(R)
print(f"diag_reg lambda*={lam_star:.6g}", file=sys.stderr)

out = pd.DataFrame({
    "gene": genes,
    "tau_unweighted": weighted_tau(mat, np.ones(len(tissues))),
    "tau_weighted_gls": censor_floor(weighted_tau(mat, w_unreg)),
    "tau_weighted_diag_reg_inv_corr": weighted_tau(mat, w_reg),
    "n_negative_weights_unregularized": n_negative,
})
out.to_csv(OUT_TAU_CSV, index=False)

weights_out = pd.DataFrame({
    "tissue": tissues,
    "weight_gls": w_unreg,
    "weight_diag_reg_inv_corr": w_reg,
    "is_negative_gls": neg_mask,
    "diag_reg_lambda": lam_star,
})
weights_out.to_csv(OUT_WEIGHTS_CSV, index=False)

print(f"genes={len(out)} tissues={len(tissues)} diag_reg_lambda={lam_star:.6g} "
      f"n_negative_unregularized={n_negative}", file=sys.stderr)
