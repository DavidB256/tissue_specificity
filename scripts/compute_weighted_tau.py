"""Weighted tau: extends tau (Yanai et al. 2005) with GLS weights that account for
cross-tissue correlation (e.g. redundant brain subregions), so correlated tissues
don't each pull full independent weight in the specificity average. R = tissue x
tissue Pearson correlation of log10(TPM+1), computed across genes. GLS weights w
solve R w = 1 (proportional to R^-1 1); w=1 for every tissue recovers unweighted tau
exactly when R=I, so no extra rescaling is needed. Negative entries in w are possible
under multicollinearity and are reported, not clipped. Ridge-shrunk weights use
R_lambda = (1-lambda)*R + lambda*I with the smallest lambda that makes every weight
positive (see gls_tau.ridge_lambda_search). Ties for a gene's max-expressing tissue
go to the first occurrence (np.argmax default)."""
import sys
import numpy as np
import pandas as pd
from gls_tau import tissue_corr, gls_weights, ridge_lambda_search, weighted_tau

IN_GCT = sys.argv[1]
OUT_TAU_CSV = sys.argv[2]
OUT_WEIGHTS_CSV = sys.argv[3]

df = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(df.columns[2:])
genes = df["Name"].values
mat = df[tissues].values.astype(float)  # genes x tissues

R = tissue_corr(mat)
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
    "tau_weighted_gls": weighted_tau(mat, w_unreg),
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
