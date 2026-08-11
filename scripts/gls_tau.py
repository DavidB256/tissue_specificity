"""Shared GLS/ridge tissue-weighting + weighted-tau helpers, used by
compute_weighted_tau.py and compute_tau_stability.py."""
import numpy as np


def tissue_corr(mat):
    """Pearson correlation of log10(TPM+1) across tissues (columns) of a genes x
    tissues matrix."""
    return np.corrcoef(np.log10(mat + 1.0), rowvar=False)


def gls_weights(R, lam=0.0):
    n = R.shape[0]
    R_lam = (1 - lam) * R + lam * np.eye(n)
    return np.linalg.solve(R_lam, np.ones(n))


def ridge_lambda_search(R, n_iter=60):
    """Smallest lambda in [0,1] s.t. gls_weights(R, lambda) is all positive.
    lambda=1 (R_lambda=I) is always all-positive, so the search is well-posed."""
    w0 = gls_weights(R, 0.0)
    if np.all(w0 > 0):
        return 0.0, w0
    lo, hi = 0.0, 1.0
    for _ in range(n_iter):
        mid = (lo + hi) / 2
        if np.all(gls_weights(R, mid) > 0):
            hi = mid
        else:
            lo = mid
    return hi, gls_weights(R, hi)


def weighted_tau(mat, w):
    """mat: genes x tissues linear TPM. w: tissue weight vector. Ties for a gene's
    max-expressing tissue go to the first occurrence (np.argmax default)."""
    max_i = mat.max(axis=1)
    argmax_idx = mat.argmax(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        term = 1 - mat / max_i[:, None]
    numer = term @ w
    denom = w.sum() - w[argmax_idx]
    return np.where(max_i > 0, numer / denom, np.nan)
