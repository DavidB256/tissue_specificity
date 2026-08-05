"""Scatterplots: gene mean-TPM vs tau, tau vs beta_posterior_median, mean-TPM vs
beta_posterior_median. OLS trendline (red) in each panel, fit in the panel's own
display scale (log10 where the axis is log-scaled) so it renders straight.
beta_posterior_median is bounded in [0,1] (a posterior probability), so its axes
are linear, unlike the old unbounded exp_beta_ratio."""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IN_TAU = sys.argv[1]
IN_POSTERIOR = sys.argv[2]
OUT_PNG = sys.argv[3]

tau_df = pd.read_csv(IN_TAU)
post_df = pd.read_csv(IN_POSTERIOR)[["gene", "beta_posterior_median"]]
df = tau_df.merge(post_df, on="gene", how="inner")
df = df.dropna(subset=["tau", "beta_posterior_median", "mean_tpm"])
df = df[df["mean_tpm"] > 0]


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
scatter_with_trend(axes[1], df["tau"], df["beta_posterior_median"], xlog=False, ylog=False,
                    xlabel="tau", ylabel="beta_posterior_median")
scatter_with_trend(axes[2], df["mean_tpm"], df["beta_posterior_median"], xlog=True, ylog=False,
                    xlabel="mean TPM across tissues (log scale)", ylabel="beta_posterior_median")

fig.suptitle(f"Gene specificity (tau) vs exp_beta P(beta|x) posterior (median across tissues, n={len(df)} genes)")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)
