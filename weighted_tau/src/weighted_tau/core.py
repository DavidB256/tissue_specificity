"""tau_weighted: tissue-specificity index tau (Yanai et al. 2005) with tissue weights
chosen by optimization.

Pipeline (aggregation and tau on log10(TPM+1) by default, tau_space="linear" for TPM):
  1. TPM < tpm_threshold -> 0; genes with no tissue >= threshold are excluded.
  2. Tissues are combined into groups (user labels, or hierarchical clustering) and
     aggregated per group (mean/max/median) into a reduced matrix. Unweighted tau on
     the reduced matrix is the reference.
  3. Per-tissue weights w on the full matrix are optimized to maximize the Spearman
     rank correlation between weighted tau (full matrix) and the reference. Optionally
     tissues in the same hierarchical cluster (share_threshold) share one weight,
     reducing the number of free parameters.
  4. Step 3 is repeated n_replicates times (seeds seed, seed+1, ...); best rho wins.

weighted tau for gene g: sum_t w_t * (1 - x_gt / max_t x_gt) / (sum_t w_t - w_argmax(g)).
Ties for the max-expressing tissue go to the first occurrence (np.argmax)."""
from __future__ import annotations

import sys
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.optimize import differential_evolution
from scipy.spatial.distance import squareform
from scipy.stats import rankdata

_AGG = {"mean": np.mean, "max": np.max, "median": np.median}


@dataclass
class TauWeightedResult:
    tau: pd.Series                  # weighted tau per gene (NaN for excluded genes)
    weights: pd.Series              # best per-tissue weights (unnormalized)
    rho: float                      # best Spearman rho vs reference
    best_replicate: int             # index into replicates
    replicates: pd.DataFrame        # per-replicate: seed, rho, n_iterations
    replicate_weights: pd.DataFrame  # tissues x replicates
    reference_tau: pd.Series        # unweighted tau on the group-aggregated matrix
    groups: pd.Series               # tissue -> group label used for the reference
    weight_clusters: pd.Series | None  # tissue -> cluster id sharing a weight (None if free)
    cluster_weights: pd.Series | None  # best weight per weight cluster (None if free)
    n_genes_used: int               # genes entering the Spearman objective


def _tissue_corr(log_mat: np.ndarray, corr_method: str) -> np.ndarray:
    if corr_method == "pearson":
        return np.corrcoef(log_mat, rowvar=False)
    if corr_method == "spearman":
        from scipy.stats import spearmanr
        return spearmanr(log_mat, axis=0)[0]
    raise ValueError(f"corr_method must be 'pearson' or 'spearman', got {corr_method!r}")


def _cluster_labels(corr: np.ndarray, threshold: float) -> np.ndarray:
    """Average linkage on 1 - r, flat cut at distance threshold (fcluster: <=)."""
    Z = linkage(squareform(1.0 - corr, checks=False), method="average")
    return fcluster(Z, t=threshold, criterion="distance")


def _tau(x: np.ndarray, w: np.ndarray) -> np.ndarray:
    max_i = x.max(axis=1)
    argmax_idx = x.argmax(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        term = 1.0 - x / max_i[:, None]
        tau = (term @ w) / (w.sum() - w[argmax_idx])
    return np.where(max_i > 0, tau, np.nan)


class _NegSpearman:
    """Picklable objective (top-level class so DE workers>1 can use it)."""

    def __init__(self, x, tau_ref, cluster_of_tissue):
        self.term = 1.0 - x / x.max(axis=1)[:, None]
        self.argmax_idx = x.argmax(axis=1)
        rank_ref = rankdata(tau_ref)
        self.rank_ref_c = rank_ref - rank_ref.mean()
        self.denom_ref = np.sqrt((self.rank_ref_c ** 2).sum())
        self.cluster_of_tissue = cluster_of_tissue

    def __call__(self, p):
        w = p if self.cluster_of_tissue is None else p[self.cluster_of_tissue]
        tau = self.term @ w / (w.sum() - w[self.argmax_idx])
        r = rankdata(tau)
        r = r - r.mean()
        return -(r @ self.rank_ref_c) / (np.sqrt((r ** 2).sum()) * self.denom_ref)


def tau_weighted(
    expr: pd.DataFrame,
    groups: Sequence | Mapping | None = None,
    group_threshold: float = 0.1,
    share_threshold: float | None = 0.1,
    corr_method: str = "pearson",
    aggregation: str = "median",
    tpm_threshold: float = 1.0,
    tau_space: str = "log",
    method: Callable = differential_evolution,
    method_kwargs: dict | None = None,
    min_weight: float = 0.01,
    max_weight: float = 1.0,
    tol: float = 1e-8,
    maxiter: int = 300,
    popsize: int = 15,
    n_replicates: int = 16,
    n_jobs: int = 8,
    seed: int = 0,
    verbose: bool = False,
) -> TauWeightedResult:
    """Compute tau_weighted for a genes x tissues TPM DataFrame.

    expr: rows = genes (or transcripts), columns = tissues, linear TPM.
    groups: how to combine tissues for the reference. A sequence of labels (one per
        column, in column order), a mapping tissue -> label (must cover all columns),
        or None (default) for hierarchical clustering of tissues (average linkage on
        1 - r, flat cut where 1 - r <= group_threshold).
    share_threshold: tissues in the same hierarchical cluster (same cut rule, on the
        full-resolution tissues) share one free weight. None = one free weight per tissue.
    corr_method: "pearson" or "spearman"; r computed on log10(TPM+1) over expressed genes
        (>= tpm_threshold in at least one tissue); used by both clusterings.
    aggregation: "mean", "max", or "median" across tissues in a group (log space).
    tpm_threshold: TPM below this -> 0; genes zero in all tissues are excluded.
    tau_space: "log" (default) aggregates and computes tau on log10(TPM+1); "linear" uses
        TPM as-is. TPM-threshold filtering and tissue-correlation clustering always use
        raw TPM / log10(TPM+1) respectively, independent of this setting.
    method: any callable with scipy.optimize.differential_evolution's signature
        method(func, bounds, seed=, tol=, **kwargs) returning an object with .x/.fun.
        For differential_evolution itself, maxiter/popsize/n_jobs (workers) apply and
        polish is off (objective is piecewise-constant, so gradient polishing is useless);
        other callables get method_kwargs only, and n_jobs is ignored.
    min_weight, max_weight: closed bounds for each weight.
    tol: convergence threshold passed as `tol`.
    n_replicates: independent optimization runs (seeds seed .. seed+n_replicates-1);
        run sequentially, each parallelized across n_jobs processes by SciPy DE
        (workers=n_jobs). Callers on Windows/macOS need an `if __name__ == "__main__":`
        guard when n_jobs > 1.
    Returns TauWeightedResult for the replicate with the highest rho."""
    if aggregation not in _AGG:
        raise ValueError(f"aggregation must be one of {sorted(_AGG)}, got {aggregation!r}")
    tissues = list(expr.columns)
    n_tissues = len(tissues)
    raw = expr.to_numpy(dtype=float)
    floored = np.where(raw >= tpm_threshold, raw, 0.0)
    expressed = floored.max(axis=1) > 0
    if tau_space not in ("log", "linear"):
        raise ValueError(f"tau_space must be 'log' or 'linear', got {tau_space!r}")
    x = np.log10(floored + 1.0) if tau_space == "log" else floored

    corr = None
    if groups is None or share_threshold is not None:
        corr = _tissue_corr(np.log10(floored[expressed] + 1.0), corr_method)

    # --- tissue groups for the reference
    if groups is None:
        labels = _cluster_labels(corr, group_threshold)
    elif isinstance(groups, Mapping):
        missing = [t for t in tissues if t not in groups]
        if missing:
            raise ValueError(f"groups mapping missing tissues: {missing}")
        labels = np.array([groups[t] for t in tissues])
    else:
        if len(groups) != n_tissues:
            raise ValueError(f"groups has {len(groups)} labels for {n_tissues} tissues")
        labels = np.asarray(groups)
    group_ids, group_idx = np.unique(labels, return_inverse=True)
    agg = _AGG[aggregation]
    x_ref = np.column_stack([agg(x[:, group_idx == g], axis=1) for g in range(len(group_ids))])
    tau_ref_all = _tau(x_ref, np.ones(x_ref.shape[1]))

    # --- weight-sharing clusters
    if share_threshold is None:
        cluster_of_tissue, n_dims = None, n_tissues
        share_labels = None
    else:
        share_labels = _cluster_labels(corr, share_threshold)
        cluster_ids, cluster_of_tissue = np.unique(share_labels, return_inverse=True)
        n_dims = len(cluster_ids)

    valid = np.isfinite(tau_ref_all) & (x.max(axis=1) > 0)
    objective = _NegSpearman(x[valid], tau_ref_all[valid], cluster_of_tissue)
    bounds = [(min_weight, max_weight)] * n_dims
    is_de = method is differential_evolution
    kwargs = dict(method_kwargs or {})
    if is_de:
        kwargs.setdefault("maxiter", maxiter)
        kwargs.setdefault("popsize", popsize)
        kwargs.setdefault("polish", False)
        kwargs.setdefault("workers", n_jobs)
        if kwargs["workers"] != 1:
            kwargs.setdefault("updating", "deferred")
    if verbose:
        print(f"genes={len(expr)} used={int(valid.sum())} tissues={n_tissues} "
              f"groups={len(group_ids)} dims={n_dims}", file=sys.stderr)

    runs = []
    for i in range(n_replicates):
        rep_kwargs = dict(kwargs)
        if verbose and is_de:
            def _cb(xk, convergence=None, _i=i):
                print(f"rep {_i} rho={-objective(xk):.6f} convergence={convergence:.4g}",
                      file=sys.stderr)
            rep_kwargs["callback"] = _cb
        res = method(objective, bounds, seed=seed + i, tol=tol, **rep_kwargs)
        runs.append((seed + i, -float(res.fun), int(getattr(res, "nit", -1)), np.asarray(res.x)))
        if verbose:
            print(f"rep {i} seed={seed + i} rho={runs[-1][1]:.6f}", file=sys.stderr)

    def to_tissue_w(p):
        return p if cluster_of_tissue is None else p[cluster_of_tissue]

    rhos = [r[1] for r in runs]
    best = int(np.argmax(rhos))
    w_best = to_tissue_w(runs[best][3])
    genes = expr.index
    return TauWeightedResult(
        tau=pd.Series(_tau(x, w_best), index=genes, name="tau_weighted"),
        weights=pd.Series(w_best, index=tissues, name="weight"),
        rho=rhos[best],
        best_replicate=best,
        replicates=pd.DataFrame({"seed": [r[0] for r in runs], "rho": rhos,
                                 "n_iterations": [r[2] for r in runs]}),
        replicate_weights=pd.DataFrame({r[0]: to_tissue_w(r[3]) for r in runs}, index=tissues),
        reference_tau=pd.Series(tau_ref_all, index=genes, name="tau_reference"),
        groups=pd.Series(labels, index=tissues, name="group"),
        weight_clusters=None if share_labels is None else pd.Series(share_labels, index=tissues, name="cluster"),
        cluster_weights=None if share_threshold is None else pd.Series(runs[best][3], index=cluster_ids, name="weight"),
        n_genes_used=int(valid.sum()),
    )
