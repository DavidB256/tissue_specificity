"""Plot exp_skewgauss2 fit for each salmon error-simulation condition
(control vs miscall_0.01/0.1/1) in a 2x2 grid."""
import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mixture_core import logpdf as comp_logpdf, unflatten_fit
from model_registry import MODELS

IN_PARQUET = sys.argv[1]
IN_FITS = sys.argv[2]
OUT_PNG = sys.argv[3]
MODEL = "exp_skewgauss2"
COMP_TYPES = MODELS[MODEL]
CONDITIONS = ["control", "miscall_0.01", "miscall_0.1", "miscall_1"]

long_df = pd.read_parquet(IN_PARQUET)
fits = pd.read_csv(IN_FITS).set_index("tissue")

fig, axes = plt.subplots(2, 2, figsize=(12, 10))
for ax, cond in zip(axes.ravel(), CONDITIONS):
    x = long_df.loc[long_df.tissue == cond, "log10_tpm"].values
    row = fits.loc[cond]
    weights, params = unflatten_fit(COMP_TYPES, row)
    ax.hist(x, bins=40, density=True, alpha=0.5, color="steelblue")
    xs = np.linspace(x.min(), x.max(), 500)
    total = np.zeros_like(xs)
    for ct, w, p in zip(COMP_TYPES, weights, params):
        comp = w * np.exp(comp_logpdf(xs, ct, p))
        total += comp
        ax.plot(xs, comp, "--", lw=1)
    ax.plot(xs, total, color="darkred", lw=2)
    ax.set_title(f"{cond}\nn={len(x)} bic={row.bic:.0f} ks={row.ks_stat:.4f} p={row.ks_pval:.2g}", fontsize=10)
    ax.set_xlabel("log10(TPM)")
fig.suptitle("exp_skewgauss2 fit vs simulated miscall rate (Kidney, SRR810007)", fontsize=13)
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)
