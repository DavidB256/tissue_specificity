"""Weighted tau using SVD/eigendecomposition tissue weights: eigendecompose R (tissue
x tissue Pearson correlation of log10(TPM+1)) as R = U diag(s) U.T, project phi = U.T @
ones(n), solve w = U @ (phi / max(s, s_floor)), normalize w /= w.sum(). s_floor=0
recovers the unregularized GLS solution R^-1 1 exactly (see gls_tau.gls_weights);
flooring small eigenvalues suppresses the blow-up they'd otherwise cause without
uniformly shrinking all of R toward I the way ridge does (see gls_tau.svd_weights).
Regularized weights use the smallest s_floor that makes every weight positive
(gls_tau.svd_floor_search, same bisection idea as compute_weighted_tau.py's
ridge_lambda_search). Contrasted with the free/free_cluster methods
(optimize_free_weights.py): svd weights are closed-form from R alone -- no numerical
optimization, no reference tau needed. TPM<1 floored to 0 before any computation. R is
restricted to genes with TPM>=1 in at least one tissue (see
gls_tau.tissue_corr_expressed_only). tau_svd (unregularized, unbounded below) is
left-censored at -1 to drop runaway outliers we don't care about."""
import sys
import numpy as np
import pandas as pd
from gls_tau import tissue_corr_expressed_only, svd_weights, svd_floor_search, weighted_tau, floor_tpm, censor_floor

IN_GCT = sys.argv[1]
OUT_TAU_CSV = sys.argv[2]
OUT_WEIGHTS_CSV = sys.argv[3]

df = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(df.columns[2:])
genes = df["Name"].values
mat = df[tissues].values.astype(float)  # genes x tissues
mat = floor_tpm(mat)  # TPM<1 -> 0

R = tissue_corr_expressed_only(mat)
w_unreg = svd_weights(R, 0.0)
neg_mask = w_unreg < 0
n_negative = int(neg_mask.sum())
if n_negative:
    print(f"negative unregularized SVD weights ({n_negative}/{len(tissues)}): "
          f"{[t for t, m in zip(tissues, neg_mask) if m]}", file=sys.stderr)

s_floor_star, w_floored = svd_floor_search(R)
print(f"s_floor*={s_floor_star:.6g}", file=sys.stderr)

out = pd.DataFrame({
    "gene": genes,
    "tau_unweighted": weighted_tau(mat, np.ones(len(tissues))),
    "tau_svd": censor_floor(weighted_tau(mat, w_unreg)),
    "tau_svd_floored": weighted_tau(mat, w_floored),
    "n_negative_weights_unregularized": n_negative,
})
out.to_csv(OUT_TAU_CSV, index=False)

weights_out = pd.DataFrame({
    "tissue": tissues,
    "weight_svd": w_unreg,
    "weight_svd_floored": w_floored,
    "is_negative_svd": neg_mask,
    "s_floor": s_floor_star,
})
weights_out.to_csv(OUT_WEIGHTS_CSV, index=False)

print(f"genes={len(out)} tissues={len(tissues)} s_floor={s_floor_star:.6g} "
      f"n_negative_unregularized={n_negative}", file=sys.stderr)
