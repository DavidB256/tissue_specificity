"""3x3 grid of fits at BIC quantiles 0,12.5,...,100% for one model."""
import sys
import json
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mixture_core import logpdf as comp_logpdf
from model_registry import MODELS

MODEL = sys.argv[1]
IN_PARQUET = sys.argv[2]
IN_FITS = sys.argv[3]
OUT_PNG = sys.argv[4]
COMP_TYPES = MODELS[MODEL]

long_df = pd.read_parquet(IN_PARQUET)
fits = pd.read_csv(IN_FITS).sort_values("bic").reset_index(drop=True)
n = len(fits)
quantiles = np.linspace(0, 1, 9)
idxs = [round(q * (n - 1)) for q in quantiles]

fig, axes = plt.subplots(3, 3, figsize=(15, 13))
for ax, q, i in zip(axes.ravel(), quantiles, idxs):
    row = fits.iloc[i]
    tissue = row.tissue
    x = long_df.loc[long_df.tissue == tissue, "log10_tpm"].values
    weights = json.loads(row.weights)
    params = json.loads(row.params)
    ax.hist(x, bins=40, density=True, alpha=0.5, color="steelblue")
    xs = np.linspace(x.min(), x.max(), 500)
    total = np.zeros_like(xs)
    for ct, w, p in zip(COMP_TYPES, weights, params):
        comp = w * np.exp(comp_logpdf(xs, ct, p))
        total += comp
        ax.plot(xs, comp, "--", lw=1)
    ax.plot(xs, total, color="darkred", lw=2)
    ax.set_title(f"BIC {q * 100:.1f}%: {tissue}\nbic={row.bic:.0f} ks={row.ks_stat:.4f}", fontsize=8)
    ax.set_xlabel("log10(TPM)", fontsize=8)
fig.suptitle(f"model={MODEL}", fontsize=13)
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)
