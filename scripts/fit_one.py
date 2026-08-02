"""Fit one named mixture model (model_registry.MODELS) to one tissue's log10_tpm. One-row CSV out."""
import sys
import json
import pandas as pd
from scipy.stats import ks_1samp
from mixture_core import fit_mixture, mixture_cdf
from model_registry import MODELS

MODEL = sys.argv[1]
TISSUE = sys.argv[2]
IN_PARQUET = sys.argv[3]
OUT_CSV = sys.argv[4]
COMP_TYPES = MODELS[MODEL]

long_df = pd.read_parquet(IN_PARQUET)
x = long_df.loc[long_df.tissue == TISSUE, "log10_tpm"].values

fit = fit_mixture(x, COMP_TYPES, n_init=5)
ks_stat, ks_p = ks_1samp(x, lambda v: mixture_cdf(v, COMP_TYPES, fit["weights"], fit["params"]))
row = {
    "tissue": TISSUE, "n": len(x), "n_params": fit["n_params"],
    "aic": fit["aic"], "bic": fit["bic"], "loglik": fit["loglik"],
    "ks_stat": ks_stat, "ks_pval": ks_p,
    "weights": json.dumps(fit["weights"]),
    "params": json.dumps(fit["params"]),
}
pd.DataFrame([row]).to_csv(OUT_CSV, index=False)
