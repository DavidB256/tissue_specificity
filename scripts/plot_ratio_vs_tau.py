"""Scatter: exp_beta_ratio_median (log10 x-axis) vs tau (y-axis), with LOESS and
quadratic-regression trend lines."""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from statsmodels.nonparametric.smoothers_lowess import lowess

IN_TAU = sys.argv[1]
IN_RATIO = sys.argv[2]
OUT_PNG = sys.argv[3]

tau_df = pd.read_csv(IN_TAU)
ratio_df = pd.read_csv(IN_RATIO)[["gene", "exp_beta_ratio_median"]]
df = tau_df.merge(ratio_df, on="gene", how="inner")
df = df.dropna(subset=["tau", "exp_beta_ratio_median", "mean_tpm"])
df = df[(df["mean_tpm"] > 0) & (df["exp_beta_ratio_median"] > 0)]

x_log = np.log10(df["exp_beta_ratio_median"].values)
y = df["tau"].values

fig, ax = plt.subplots(figsize=(8, 6))
ax.scatter(df["exp_beta_ratio_median"], y, s=2, alpha=0.2)
ax.set_xscale("log")
ax.set_xlabel("exp_beta_ratio_median (log scale)")
ax.set_ylabel("tau")

# quadratic regression (fit in log10(x) space, plotted back on the log-scaled axis)
coeffs = np.polyfit(x_log, y, 2)
xs_log = np.linspace(x_log.min(), x_log.max(), 300)
ys_quad = np.polyval(coeffs, xs_log)
ax.plot(10 ** xs_log, ys_quad, color="red", lw=2, label="quadratic fit")

# LOESS
smoothed = lowess(y, x_log, frac=0.3, return_sorted=True)
ax.plot(10 ** smoothed[:, 0], smoothed[:, 1], color="darkorange", lw=2, label="LOESS")

ax.legend()
fig.suptitle(f"tau vs exp_beta_ratio_median (n={len(df)} genes)")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)
