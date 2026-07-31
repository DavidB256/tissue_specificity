"""Fit k=1,2,3-component Gaussian mixtures per tissue on log10_tpm. Output fit table."""
import sys
import pandas as pd
import numpy as np
from sklearn.mixture import GaussianMixture
from scipy.stats import ks_1samp, norm

IN_PARQUET = sys.argv[1]
OUT_CSV = sys.argv[2]
KS = [1, 2, 3]

long_df = pd.read_parquet(IN_PARQUET)

def mixture_cdf(x, means, covs, weights):
    x = np.atleast_1d(x)
    cdf = np.zeros_like(x, dtype=float)
    for m, c, w in zip(means, covs, weights):
        cdf += w * norm.cdf(x, loc=m, scale=np.sqrt(c))
    return cdf

rows = []
for tissue, sub in long_df.groupby("tissue"):
    x = sub["log10_tpm"].values.reshape(-1, 1)
    for k in KS:
        gmm = GaussianMixture(n_components=k, random_state=0, n_init=3).fit(x)
        means = gmm.means_.ravel()
        covs = gmm.covariances_.ravel()
        weights = gmm.weights_.ravel()
        loglik = gmm.score(x) * len(x)
        ks_stat, ks_p = ks_1samp(x.ravel(), lambda v: mixture_cdf(v, means, covs, weights))
        rows.append({
            "tissue": tissue, "k": k, "n": len(x),
            "aic": gmm.aic(x), "bic": gmm.bic(x), "loglik": loglik,
            "ks_stat": ks_stat, "ks_pval": ks_p,
            "means": ";".join(f"{m:.4f}" for m in means),
            "stds": ";".join(f"{np.sqrt(c):.4f}" for c in covs),
            "weights": ";".join(f"{w:.4f}" for w in weights),
        })

out = pd.DataFrame(rows)
out.to_csv(OUT_CSV, index=False)
print(out.groupby("k")[["aic", "bic", "ks_stat"]].mean(), file=sys.stderr)
