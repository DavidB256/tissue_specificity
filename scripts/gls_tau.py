"""Shared diag_reg_inv_corr ("diagonal-regularized inverse correlation") tissue-
weighting + weighted-tau helpers, used by compute_weighted_tau.py and
compute_tau_stability.py."""
import numpy as np


def tissue_corr(mat):
    """Pearson correlation of log10(TPM+1) across tissues (columns) of a genes x
    tissues matrix."""
    return np.corrcoef(np.log10(mat + 1.0), rowvar=False)


def tissue_corr_pairwise_masked(mat):
    """Pairwise-masked variant of tissue_corr: each tissue-pair's correlation (log10
    TPM+1) is restricted to genes with TPM>=1 (post floor_tpm) in at least one of the
    two tissues, rather than one fixed gene set for the whole matrix. Violates the
    single-consistent-sample-matrix assumption behind PSD -- R here is NOT guaranteed
    PSD; downstream diag_reg_inv_corr_weights/diag_reg_lambda_search still work
    (diag_reg_lambda_search pushes lambda up until R_lam is invertible with
    all-positive weights)."""
    n = mat.shape[1]
    log_mat = np.log10(mat + 1.0)
    R = np.eye(n)
    for i in range(n):
        for j in range(i + 1, n):
            mask = (mat[:, i] > 0) | (mat[:, j] > 0)
            r = np.corrcoef(log_mat[mask, i], log_mat[mask, j])[0, 1]
            R[i, j] = R[j, i] = r
    return R


def tissue_corr_expressed_only(mat):
    """tissue_corr restricted to genes with TPM>=1 (post floor_tpm) in at least one
    tissue. Single consistent gene set across all tissue pairs, so R stays PSD (unlike
    tissue_corr_pairwise_masked)."""
    expressed = mat.max(axis=1) > 0
    return tissue_corr(mat[expressed])


def diag_reg_inv_corr_weights(R, lam=0.0):
    """Solve R_lambda w = 1, R_lambda = (1-lam)*R + lam*I. lam=0 is the plain inverse-
    correlation solve (R^-1 1) -- routinely has negative entries under multicollinearity
    and is only useful here as the starting point diag_reg_lambda_search searches up
    from, not as a result in its own right. Raises np.linalg.LinAlgError if R_lam is
    exactly singular (e.g. lam=0 with duplicate tissue columns from resampling with
    replacement)."""
    n = R.shape[0]
    R_lam = (1 - lam) * R + lam * np.eye(n)
    return np.linalg.solve(R_lam, np.ones(n))


def diag_reg_lambda_search(R, n_iter=60):
    """Smallest lambda in [0,1] s.t. diag_reg_inv_corr_weights(R, lambda) is all
    positive -- the diagonal-regularization strength for diag_reg_inv_corr weighting
    (R_lambda = (1-lambda)*R + lambda*I, see diag_reg_inv_corr_weights). lambda=1
    (R_lambda=I) is always all-positive, so the search is well-posed. Singular R_lam
    (e.g. from duplicate tissues under resampling) is treated as "not all positive" so
    the search pushes lambda up until it's invertible."""
    def all_positive(lam):
        try:
            return np.all(diag_reg_inv_corr_weights(R, lam) > 0)
        except np.linalg.LinAlgError:
            return False

    if all_positive(0.0):
        return 0.0, diag_reg_inv_corr_weights(R, 0.0)
    lo, hi = 0.0, 1.0
    for _ in range(n_iter):
        mid = (lo + hi) / 2
        if all_positive(mid):
            hi = mid
        else:
            lo = mid
    return hi, diag_reg_inv_corr_weights(R, hi)


def svd_weights(R, s_floor=0.0):
    """SVD/eigendecomposition weights: w = U @ (phi / max(s, s_floor)), phi = U.T @ 1,
    R = U diag(s) U.T. Equivalent to diag_reg_inv_corr_weights(R, 0) when s_floor=0 (both solve R w =
    1 via the eigenbasis vs. a linear solve) but lets small eigenvalues (tight/redundant
    tissue clusters) be floored independently rather than uniformly shrinking all of R
    toward I like diagonal regularization does. s_floor=s.max() collapses to uniform weights (every
    eigenvalue floored to the same value, and U @ U.T @ 1 = 1 since U is a complete
    orthonormal basis)."""
    n = R.shape[0]
    U, s, _ = np.linalg.svd(R)
    phi = U.T @ np.ones(n)
    w = U @ (phi / np.maximum(s, s_floor))
    return w / w.sum()


def svd_floor_search(R, n_iter=60):
    """Smallest s_floor in [0, s.max()] s.t. svd_weights(R, s_floor) is all positive.
    s_floor=s.max() is always all-positive (uniform weights), so the search is
    well-posed -- same bisection structure as diag_reg_lambda_search."""
    n = R.shape[0]
    U, s, _ = np.linalg.svd(R)
    smax = s.max()

    def weights_for(floor):
        w = U @ (phi_for(floor))
        return w / w.sum()

    phi = U.T @ np.ones(n)

    def phi_for(floor):
        return phi / np.maximum(s, floor)

    def all_positive(floor):
        return np.all(weights_for(floor) > 0)

    if all_positive(0.0):
        return 0.0, weights_for(0.0)
    lo, hi = 0.0, smax
    for _ in range(n_iter):
        mid = (lo + hi) / 2
        if all_positive(mid):
            hi = mid
        else:
            lo = mid
    return hi, weights_for(hi)


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


def margin_weights(mat):
    """Margin-weighted exclusivity: per gene, margin = (max - second_max) / max,
    credited entirely to that gene's max-expressing tissue (same argmax convention as
    weighted_tau -- ties go to the first occurrence). w_t = sum over genes with
    argmax==t of margin(g), normalized to sum to 1. A tissue redundant with others
    (e.g. one of several correlated brain subregions) rarely wins a gene by a wide
    margin -- a correlated partner is usually close behind -- so it accumulates little
    weight; a distinctive tissue (e.g. testis) wins its genes decisively and
    accumulates more. Unlike diag_reg_inv_corr_weights, this never touches the tissue
    x tissue correlation matrix at all -- it's a summary of the expression matrix
    itself -- so there is no matrix inversion, no regularization search, and no
    negative-weight failure mode to guard against: margin(g) in [0,1] by construction,
    so w is elementwise >= 0 for free. Genes with max<=0 (unexpressed everywhere)
    contribute 0 to every tissue."""
    n_genes, n_tissues = mat.shape
    top2 = -np.sort(-mat, axis=1)[:, :2]
    max1 = top2[:, 0]
    max2 = top2[:, 1] if n_tissues > 1 else np.zeros(n_genes)
    with np.errstate(invalid="ignore", divide="ignore"):
        margin = np.where(max1 > 0, (max1 - max2) / max1, 0.0)
    argmax_idx = mat.argmax(axis=1)
    w = np.zeros(n_tissues)
    np.add.at(w, argmax_idx, margin)
    return w / w.sum()


def floor_tpm(mat, floor=1.0):
    """TPM<floor -> 0 (treated as unexpressed noise), matches TPM>=1 convention used
    elsewhere in the pipeline (e.g. build_expr_matrix_allgenes.py)."""
    return np.where(mat >= floor, mat, 0.0)


def censor_floor(tau, floor=-1.0):
    """Left-censor: tau < floor -> floor. Ignores NaN. Used for unregularized/unbounded-
    below weighted-tau variants (e.g. compute_svd_weighted_tau.py's tau_svd), which can
    produce runaway outliers."""
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
