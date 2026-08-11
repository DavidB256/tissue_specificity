"""Weighted tau: extends tau (Yanai et al. 2005) with GLS weights that account for
cross-tissue correlation (e.g. redundant brain subregions), so correlated tissues
don't each pull full independent weight in the specificity average. R = tissue x
tissue Pearson correlation of log10(TPM+1), computed across genes. GLS weights w
solve R w = 1 (proportional to R^-1 1); w=1 for every tissue recovers unweighted tau
exactly when R=I, so no extra rescaling is needed. Negative entries in w are possible
under multicollinearity and are reported, not clipped. Ridge-shrunk weights use
R_lambda = (1-lambda)*R + lambda*I with the smallest lambda (binary search over
[0,1]) that makes every weight positive -- lambda=1 (R_lambda=I) is always positive,
so the search is well-posed. Ties for a gene's max-expressing tissue go to the first
occurrence (np.argmax default)."""
import sys
import numpy as np
import pandas as pd

IN_GCT = sys.argv[1]
OUT_TAU_CSV = sys.argv[2]
OUT_WEIGHTS_CSV = sys.argv[3]

df = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(df.columns[2:])
genes = df["Name"].values
mat = df[tissues].values.astype(float)  # genes x tissues
n = mat.shape[1]
ones = np.ones(n)

logmat = np.log10(mat + 1.0)
R = np.corrcoef(logmat, rowvar=False)


def gls_weights(lam):
    R_lam = (1 - lam) * R + lam * np.eye(n)
    return np.linalg.solve(R_lam, ones)


w_unreg = gls_weights(0.0)
neg_mask = w_unreg < 0
n_negative = int(neg_mask.sum())
if n_negative:
    print(f"negative unregularized GLS weights ({n_negative}/{n}): "
          f"{[t for t, m in zip(tissues, neg_mask) if m]}", file=sys.stderr)

# binary search: smallest lambda in [0,1] s.t. all(gls_weights(lambda) > 0)
w_hi = gls_weights(1.0)
assert np.all(w_hi > 0), "ridge weights at lambda=1 (identity) must be all positive"
if n_negative == 0:
    lam_star, w_reg = 0.0, w_unreg
else:
    lo, hi = 0.0, 1.0
    for _ in range(60):
        mid = (lo + hi) / 2
        if np.all(gls_weights(mid) > 0):
            hi = mid
        else:
            lo = mid
    lam_star = hi
    w_reg = gls_weights(lam_star)
print(f"ridge lambda*={lam_star:.6g}", file=sys.stderr)

max_i = mat.max(axis=1)
argmax_idx = mat.argmax(axis=1)  # ties -> first occurrence
with np.errstate(invalid="ignore", divide="ignore"):
    term = 1 - mat / max_i[:, None]  # genes x tissues, 0 at each gene's max tissue


def weighted_tau(w):
    numer = term @ w
    denom = w.sum() - w[argmax_idx]  # exclude each gene's own max-tissue weight
    return np.where(max_i > 0, numer / denom, np.nan)


out = pd.DataFrame({
    "gene": genes,
    "tau_unweighted": weighted_tau(ones),
    "tau_weighted_gls": weighted_tau(w_unreg),
    "tau_weighted_gls_ridge": weighted_tau(w_reg),
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

print(f"genes={len(out)} tissues={n} ridge_lambda={lam_star:.6g} "
      f"n_negative_unregularized={n_negative}", file=sys.stderr)
