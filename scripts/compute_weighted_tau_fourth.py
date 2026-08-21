"""Weighted tau using the elementwise-4th-power tissue correlation matrix R**4 (Hadamard
power -- stays PSD via repeated Schur product theorem, diagonal stays 1, up-weights very
high |corr| much more sharply than R**2; valid since all entries of R are non-negative
here, so no sign flip). Otherwise identical to compute_weighted_tau_squared.py: GLS
weights solve R_4 w = 1, ridge-shrunk to the smallest lambda making every weight positive."""
import sys
import numpy as np
import pandas as pd
from gls_tau import tissue_corr_expressed_only, gls_weights, ridge_lambda_search, weighted_tau, floor_tpm, censor_floor

IN_GCT = sys.argv[1]
OUT_TAU_CSV = sys.argv[2]
OUT_WEIGHTS_CSV = sys.argv[3]

df = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(df.columns[2:])
genes = df["Name"].values
mat = df[tissues].values.astype(float)  # genes x tissues
mat = floor_tpm(mat)  # TPM<1 -> 0

R = tissue_corr_expressed_only(mat)
assert (R >= 0).all(), "raising to a power flips sign of negative correlations -- not valid here"
R_4 = R ** 4

w_unreg = gls_weights(R_4, 0.0)
neg_mask = w_unreg < 0
n_negative = int(neg_mask.sum())
if n_negative:
    print(f"negative unregularized GLS weights ({n_negative}/{len(tissues)}): "
          f"{[t for t, m in zip(tissues, neg_mask) if m]}", file=sys.stderr)

lam_star, w_reg = ridge_lambda_search(R_4)
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
