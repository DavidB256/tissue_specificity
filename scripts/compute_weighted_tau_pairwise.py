"""Pairwise-masked variant of compute_weighted_tau.py: R comes from
gls_tau.tissue_corr_pairwise_masked (each tissue-pair's correlation restricted to
genes with TPM>=1 in at least one of the two tissues) instead of gls_tau.tissue_corr.
R is not guaranteed PSD here; ridge_lambda_search still finds a lambda with all-positive
weights since R_lambda=I (lambda=1) is always valid, but unregularized (lambda=0)
weights lose their usual GLS interpretation. See gls_tau.tissue_corr_pairwise_masked."""
import sys
import numpy as np
import pandas as pd
from gls_tau import tissue_corr_pairwise_masked, gls_weights, ridge_lambda_search, weighted_tau, floor_tpm, censor_floor

IN_GCT = sys.argv[1]
OUT_TAU_CSV = sys.argv[2]
OUT_WEIGHTS_CSV = sys.argv[3]

df = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(df.columns[2:])
genes = df["Name"].values
mat = df[tissues].values.astype(float)  # genes x tissues
mat = floor_tpm(mat)  # TPM<1 -> 0

R = tissue_corr_pairwise_masked(mat)
w_unreg = gls_weights(R, 0.0)
neg_mask = w_unreg < 0
n_negative = int(neg_mask.sum())
if n_negative:
    print(f"negative unregularized GLS weights ({n_negative}/{len(tissues)}): "
          f"{[t for t, m in zip(tissues, neg_mask) if m]}", file=sys.stderr)

lam_star, w_reg = ridge_lambda_search(R)
print(f"ridge lambda*={lam_star:.6g}", file=sys.stderr)

out = pd.DataFrame({
    "gene": genes,
    "tau_unweighted": weighted_tau(mat, np.ones(len(tissues))),
    "tau_weighted_gls": censor_floor(weighted_tau(mat, w_unreg)),
    "tau_weighted_gls_ridge": weighted_tau(mat, w_reg),
    "n_negative_weights_unregularized": n_negative,
})
out.to_csv(OUT_TAU_CSV, index=False)

weights_out = pd.DataFrame({
    "tissue": tissues,
    "weight_gls": w_unreg,
    "weight_gls_ridge": w_reg,
    "is_negative_gls": neg_mask,
    "ridge_lambda": lam_star,
})
weights_out.to_csv(OUT_WEIGHTS_CSV, index=False)

print(f"genes={len(out)} tissues={len(tissues)} ridge_lambda={lam_star:.6g} "
      f"n_negative_unregularized={n_negative}", file=sys.stderr)
