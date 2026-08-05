"""Scatterplots: gene mean-TPM vs tau, tau vs exp_beta_ratio_median, mean-TPM vs
exp_beta_ratio_median. OLS trendline (red) in each panel, fit in the panel's own
display scale (log10 where the axis is log-scaled) so it renders straight."""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IN_TAU = sys.argv[1]
IN_RATIO = sys.argv[2]
OUT_PNG = sys.argv[3]

tau_df = pd.read_csv(IN_TAU)
ratio_df = pd.read_csv(IN_RATIO)[["gene", "exp_beta_ratio_median"]]
df = tau_df.merge(ratio_df, on="gene", how="inner")
df = df.dropna(subset=["tau", "exp_beta_ratio_median", "mean_tpm"])
df = df[(df["mean_tpm"] > 0) & (df["exp_beta_ratio_median"] > 0)]


def scatter_with_trend(ax, x, y, xlog, ylog, xlabel, ylabel):
    ax.scatter(x, y, s=2, alpha=0.2)
    if xlog:
        ax.set_xscale("log")
    if ylog:
        ax.set_yscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    xt = np.log10(x) if xlog else x
    yt = np.log10(y) if ylog else y
    slope, intercept = np.polyfit(xt, yt, 1)
    xt_line = np.linspace(xt.min(), xt.max(), 200)
    yt_line = slope * xt_line + intercept
    x_line = 10 ** xt_line if xlog else xt_line
    y_line = 10 ** yt_line if ylog else yt_line
    ax.plot(x_line, y_line, color="red", lw=2)


fig, axes = plt.subplots(1, 3, figsize=(18, 5))

scatter_with_trend(axes[0], df["mean_tpm"], df["tau"], xlog=True, ylog=False,
                    xlabel="mean TPM across tissues (log scale)", ylabel="tau")
scatter_with_trend(axes[1], df["tau"], df["exp_beta_ratio_median"], xlog=False, ylog=True,
                    xlabel="tau", ylabel="exp_beta_ratio_median (log scale)")
scatter_with_trend(axes[2], df["mean_tpm"], df["exp_beta_ratio_median"], xlog=True, ylog=True,
                    xlabel="mean TPM across tissues (log scale)", ylabel="exp_beta_ratio_median (log scale)")

fig.suptitle(f"Gene specificity (tau) vs exp_beta exponential/beta likelihood ratio (median across tissues, n={len(df)} genes)")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)
