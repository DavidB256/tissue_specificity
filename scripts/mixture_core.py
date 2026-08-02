"""Generic 1D mixture-model fitting via direct MLE (scipy.optimize), for arbitrary
combinations of norm / expon / skewnorm / t components."""
import numpy as np
from scipy.stats import norm, expon, skewnorm, t as student_t
from scipy.special import logsumexp
from scipy.optimize import minimize

COMP_NPARAMS = {"norm": 2, "expon": 1, "skewnorm": 3, "t": 3}


def n_params(comp_types):
    return sum(COMP_NPARAMS[c] for c in comp_types) + (len(comp_types) - 1)


def _raw_to_params(ct, raw):
    if ct == "expon":
        return {"rate": float(np.exp(raw[0]))}
    if ct == "norm":
        return {"mean": float(raw[0]), "std": float(np.exp(raw[1]))}
    if ct == "skewnorm":
        return {"loc": float(raw[0]), "scale": float(np.exp(raw[1])), "a": float(raw[2])}
    if ct == "t":
        return {"loc": float(raw[0]), "scale": float(np.exp(raw[1])), "df": float(np.clip(np.exp(np.clip(raw[2], -50, 50)), 0.3, 50))}
    raise ValueError(ct)


def _params_to_raw(ct, params):
    if ct == "expon":
        return [np.log(params["rate"])]
    if ct == "norm":
        return [params["mean"], np.log(params["std"])]
    if ct == "skewnorm":
        return [params["loc"], np.log(params["scale"]), params["a"]]
    if ct == "t":
        return [params["loc"], np.log(params["scale"]), np.log(params["df"])]
    raise ValueError(ct)


def _unpack(theta, comp_types):
    k = len(comp_types)
    logits = theta[: k - 1]
    w = np.exp(np.concatenate([logits, [0.0]]))
    weights = w / w.sum()
    idx = k - 1
    params = []
    for ct in comp_types:
        npar = COMP_NPARAMS[ct]
        params.append(_raw_to_params(ct, theta[idx: idx + npar]))
        idx += npar
    return weights, params


def canonicalize_by_weight(comp_types, weights, params_list):
    """Fix label-switching: within each group of identically-typed components (e.g.
    the two normals in gauss2), reorder by descending weight. Component order among
    same-type components is otherwise arbitrary (an artifact of which random restart
    the optimizer landed on), so without this, two fits with essentially the same
    density can land far apart in the raw parameter space (theta_from_fit references
    weight logits against whichever component ends up last)."""
    groups = {}
    for i, ct in enumerate(comp_types):
        groups.setdefault(ct, []).append(i)
    new_weights = list(weights)
    new_params = list(params_list)
    for idxs in groups.values():
        if len(idxs) < 2:
            continue
        order = sorted(idxs, key=lambda i: -weights[i])
        for slot, src in zip(idxs, order):
            new_weights[slot] = weights[src]
            new_params[slot] = params_list[src]
    return new_weights, new_params


def theta_from_fit(comp_types, weights, params_list):
    """Reconstruct the raw unconstrained optimization vector (k-1 weight logits,
    referenced against the last component, plus each component's raw params) from a
    stored fit's constrained weights/params. Inverse of _unpack, used to place tissues
    in the model's actual fitting parameter space (e.g. for PCA)."""
    w_last = max(weights[-1], 1e-12)
    logits = [np.log(max(w, 1e-12) / w_last) for w in weights[:-1]]
    raw = list(logits)
    for ct, p in zip(comp_types, params_list):
        raw += _params_to_raw(ct, p)
    return np.array(raw, dtype=float)


def logpdf(x, ct, params):
    if ct == "expon":
        return expon.logpdf(x, scale=1.0 / params["rate"])
    if ct == "norm":
        return norm.logpdf(x, loc=params["mean"], scale=params["std"])
    if ct == "skewnorm":
        return skewnorm.logpdf(x, a=params["a"], loc=params["loc"], scale=params["scale"])
    if ct == "t":
        return student_t.logpdf(x, df=params["df"], loc=params["loc"], scale=params["scale"])
    raise ValueError(ct)


def cdf(x, ct, params):
    if ct == "expon":
        return expon.cdf(x, scale=1.0 / params["rate"])
    if ct == "norm":
        return norm.cdf(x, loc=params["mean"], scale=params["std"])
    if ct == "skewnorm":
        return skewnorm.cdf(x, a=params["a"], loc=params["loc"], scale=params["scale"])
    if ct == "t":
        return student_t.cdf(x, df=params["df"], loc=params["loc"], scale=params["scale"])
    raise ValueError(ct)


def mixture_logpdf(x, comp_types, weights, params_list):
    logmat = np.stack([np.log(w) + logpdf(x, ct, p) for ct, w, p in zip(comp_types, weights, params_list)], axis=1)
    return logsumexp(logmat, axis=1)


def mixture_cdf(x, comp_types, weights, params_list):
    return sum(w * cdf(x, ct, p) for ct, w, p in zip(comp_types, weights, params_list))


def _neg_loglik(theta, x, comp_types):
    weights, params = _unpack(theta, comp_types)
    ll = mixture_logpdf(x, comp_types, weights, params)
    if not np.all(np.isfinite(ll)):
        return 1e12
    return -ll.sum()


def _init_theta(x, comp_types, rng):
    k = len(comp_types)
    logits = rng.normal(0, 0.5, size=k - 1)
    q = np.quantile(x, np.linspace(0.15, 0.85, k))
    raw = []
    for ct, center in zip(comp_types, q):
        center = max(center, 0.05)
        if ct == "expon":
            p = {"rate": 1.0 / center * rng.uniform(0.7, 1.3)}
        elif ct == "norm":
            p = {"mean": center * rng.uniform(0.8, 1.2), "std": max(x.std() * rng.uniform(0.3, 0.7), 0.05)}
        elif ct == "skewnorm":
            p = {"loc": center * rng.uniform(0.8, 1.2), "scale": max(x.std() * rng.uniform(0.3, 0.7), 0.05), "a": rng.uniform(-3, 3)}
        elif ct == "t":
            p = {"loc": center * rng.uniform(0.8, 1.2), "scale": max(x.std() * rng.uniform(0.3, 0.7), 0.05), "df": rng.uniform(3, 15)}
        raw += _params_to_raw(ct, p)
    return np.concatenate([logits, raw])


def fit_mixture(x, comp_types, n_init=5, seed0=0):
    x = np.asarray(x, dtype=float)
    best = None
    for i in range(n_init):
        rng = np.random.RandomState(seed0 + i)
        theta0 = _init_theta(x, comp_types, rng)
        res = minimize(_neg_loglik, theta0, args=(x, comp_types), method="L-BFGS-B")
        if not np.isfinite(res.fun):
            continue
        if best is None or res.fun < best.fun:
            best = res
    weights, params = _unpack(best.x, comp_types)
    loglik = -best.fun
    npar = n_params(comp_types)
    n = len(x)
    aic = 2 * npar - 2 * loglik
    bic = npar * np.log(n) - 2 * loglik
    return {"weights": [float(w) for w in weights], "params": params, "loglik": loglik,
            "n_params": npar, "aic": aic, "bic": bic}
