"""Generic EM for mixtures of Exponential and Gaussian components on 1D data."""
import numpy as np
from scipy.stats import expon, norm

MIN_STD = 1e-2
MIN_RATE = 1e-3
EPS = 1e-10


def _pdf(x, ctype, params):
    if ctype == "expon":
        return expon.pdf(x, scale=1.0 / params["rate"])
    return norm.pdf(x, loc=params["mean"], scale=params["std"])


def _cdf(x, ctype, params):
    if ctype == "expon":
        return expon.cdf(x, scale=1.0 / params["rate"])
    return norm.cdf(x, loc=params["mean"], scale=params["std"])


def mixture_cdf(x, ctypes, weights, params_list):
    x = np.atleast_1d(x)
    cdf = np.zeros_like(x, dtype=float)
    for ctype, w, params in zip(ctypes, weights, params_list):
        cdf += w * _cdf(x, ctype, params)
    return cdf


def _init_params(x, ctypes, rng):
    q = np.quantile(x, np.linspace(0.15, 0.85, len(ctypes)))
    params_list = []
    for ctype, center in zip(ctypes, q):
        if ctype == "expon":
            m = max(center, 0.05)
            params_list.append({"rate": 1.0 / m * rng.uniform(0.7, 1.3)})
        else:
            params_list.append({"mean": center * rng.uniform(0.8, 1.2), "std": x.std() * rng.uniform(0.3, 0.7)})
    weights = rng.dirichlet(np.ones(len(ctypes)))
    return weights, params_list


def _e_step(x, ctypes, weights, params_list):
    dens = np.stack([w * _pdf(x, ct, p) for ct, w, p in zip(ctypes, weights, params_list)], axis=1)
    total = dens.sum(axis=1, keepdims=True) + EPS
    resp = dens / total
    loglik = np.log(total.ravel() + EPS).sum()
    return resp, loglik


def _m_step(x, ctypes, resp):
    n = len(x)
    weights = resp.sum(axis=0) / n
    params_list = []
    for k, ctype in enumerate(ctypes):
        r = resp[:, k]
        rsum = r.sum() + EPS
        if ctype == "expon":
            mean = (r * x).sum() / rsum
            rate = 1.0 / max(mean, MIN_RATE)
            params_list.append({"rate": rate})
        else:
            mean = (r * x).sum() / rsum
            var = (r * (x - mean) ** 2).sum() / rsum
            params_list.append({"mean": mean, "std": max(np.sqrt(var), MIN_STD)})
    return weights, params_list


def fit_em(x, ctypes, n_init=5, max_iter=300, tol=1e-6):
    """Fit a mixture of components (ctypes: list of 'expon'/'norm') via EM. Returns best-of-n_init fit."""
    x = np.asarray(x, dtype=float)
    best = None
    for seed in range(n_init):
        rng = np.random.RandomState(seed)
        weights, params_list = _init_params(x, ctypes, rng)
        prev_ll = -np.inf
        for _ in range(max_iter):
            resp, ll = _e_step(x, ctypes, weights, params_list)
            weights, params_list = _m_step(x, ctypes, resp)
            if abs(ll - prev_ll) < tol * abs(prev_ll if prev_ll != 0 else 1):
                prev_ll = ll
                break
            prev_ll = ll
        if best is None or prev_ll > best["loglik"]:
            best = {"weights": weights, "params_list": params_list, "loglik": prev_ll}
    n_params = sum(1 if ct == "expon" else 2 for ct in ctypes) + (len(ctypes) - 1)
    n = len(x)
    aic = 2 * n_params - 2 * best["loglik"]
    bic = n_params * np.log(n) - 2 * best["loglik"]
    return {**best, "n_params": n_params, "aic": aic, "bic": bic}
