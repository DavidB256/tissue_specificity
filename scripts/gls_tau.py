"""Shared GLS/ridge tissue-weighting + weighted-tau helpers, used by
compute_weighted_tau.py and compute_tau_stability.py."""
import numpy as np


def tissue_corr(mat):
    """Pearson correlation of log10(TPM+1) across tissues (columns) of a genes x
    tissues matrix."""
    return np.corrcoef(np.log10(mat + 1.0), rowvar=False)


def gls_weights(R, lam=0.0):
    """Raises np.linalg.LinAlgError if R_lam is exactly singular (e.g. lam=0 with
    duplicate tissue columns from resampling with replacement)."""
    n = R.shape[0]
    R_lam = (1 - lam) * R + lam * np.eye(n)
    return np.linalg.solve(R_lam, np.ones(n))


def ridge_lambda_search(R, n_iter=60):
    """Smallest lambda in [0,1] s.t. gls_weights(R, lambda) is all positive.
    lambda=1 (R_lambda=I) is always all-positive, so the search is well-posed.
    Singular R_lam (e.g. from duplicate tissues under resampling) is treated as
    "not all positive" so the search pushes lambda up until it's invertible."""
    def all_positive(lam):
        try:
            return np.all(gls_weights(R, lam) > 0)
        except np.linalg.LinAlgError:
            return False

    if all_positive(0.0):
        return 0.0, gls_weights(R, 0.0)
    lo, hi = 0.0, 1.0
    for _ in range(n_iter):
        mid = (lo + hi) / 2
        if all_positive(mid):
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


def floor_tpm(mat, floor=1.0):
    """TPM<floor -> 0 (treated as unexpressed noise), matches TPM>=1 convention used
    elsewhere in the pipeline (e.g. build_expr_matrix_allgenes.py)."""
    return np.where(mat >= floor, mat, 0.0)


def censor_floor(tau, floor=-1.0):
    """Left-censor: tau < floor -> floor. Ignores NaN. Used for tau_weighted_gls
    (unregularized), which is unbounded below and can produce runaway outliers."""
    return np.where(np.isnan(tau), tau, np.maximum(tau, floor))


GENE_SUBSETS = ["all", "highspec", "lowspec", "multitissue", "onlyneural", "maxneural"]


def gene_subset_mask(subset, mat, tissues, tau_unweighted):
    """Boolean gene mask for a named subset. mat: genes x tissues, floor_tpm'd (TPM<1 ->
    0). tissues: column names aligned to mat. tau_unweighted: per-gene tau_unweighted
    from weighted_tau.csv, aligned to mat's gene rows."""
    if subset == "all":
        return np.ones(mat.shape[0], dtype=bool)
    if subset == "highspec":
        return tau_unweighted > 0.9
    if subset == "lowspec":
        return tau_unweighted <= 0.6
    if subset == "multitissue":
        return (mat > 0).sum(axis=1) > 1
    brain_idx = np.array([i for i, t in enumerate(tissues) if "Brain" in t])
    nonbrain_idx = np.array([i for i in range(len(tissues)) if i not in brain_idx])
    if subset == "onlyneural":
        return (mat[:, nonbrain_idx] > 0).sum(axis=1) == 0
    if subset == "maxneural":
        return np.isin(mat.argmax(axis=1), brain_idx)
    raise ValueError(f"unknown gene subset: {subset}")
