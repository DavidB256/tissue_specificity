"""Fit Exponential + 2-Gaussian mixture per tissue on log10_tpm. Output fit table."""
import sys
import pandas as pd
import numpy as np
from scipy.stats import ks_1samp
from mixture_em import fit_em, mixture_cdf

IN_PARQUET = sys.argv[1]
OUT_CSV = sys.argv[2]
CTYPES = ["expon", "norm", "norm"]

long_df = pd.read_parquet(IN_PARQUET)

rows = []
for tissue, sub in long_df.groupby("tissue"):
    x = sub["log10_tpm"].values
    fit = fit_em(x, CTYPES, n_init=5)
    weights, params_list = fit["weights"], fit["params_list"]
    ks_stat, ks_p = ks_1samp(x, lambda v: mixture_cdf(v, CTYPES, weights, params_list))
    rate = params_list[0]["rate"]
    m1, s1 = params_list[1]["mean"], params_list[1]["std"]
    m2, s2 = params_list[2]["mean"], params_list[2]["std"]
    rows.append({
        "tissue": tissue, "n": len(x),
        "aic": fit["aic"], "bic": fit["bic"], "loglik": fit["loglik"],
        "ks_stat": ks_stat, "ks_pval": ks_p,
        "w_expon": weights[0], "rate": rate,
        "w_norm1": weights[1], "mean1": m1, "std1": s1,
        "w_norm2": weights[2], "mean2": m2, "std2": s2,
    })

out = pd.DataFrame(rows)
out.to_csv(OUT_CSV, index=False)
print(out[["aic", "bic", "ks_stat"]].mean(), file=sys.stderr)
