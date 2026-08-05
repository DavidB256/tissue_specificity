"""Scatterplots: gene mean-TPM vs tau, tau vs exp_beta_ratio, mean-TPM vs exp_beta_ratio."""
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IN_TAU = sys.argv[1]
IN_RATIO = sys.argv[2]
OUT_PNG = sys.argv[3]

tau_df = pd.read_csv(IN_TAU)
ratio_df = pd.read_csv(IN_RATIO)[["gene", "exp_beta_ratio"]]
df = tau_df.merge(ratio_df, on="gene", how="inner")
df = df.dropna(subset=["tau", "exp_beta_ratio", "mean_tpm"])
df = df[(df["mean_tpm"] > 0) & (df["exp_beta_ratio"] > 0)]

fig, axes = plt.subplots(1, 3, figsize=(18, 5))

axes[0].scatter(df["mean_tpm"], df["tau"], s=2, alpha=0.2)
axes[0].set_xscale("log")
axes[0].set_xlabel("mean TPM across tissues (log scale)")
axes[0].set_ylabel("tau")

axes[1].scatter(df["tau"], df["exp_beta_ratio"], s=2, alpha=0.2)
axes[1].set_yscale("log")
axes[1].set_xlabel("tau")
axes[1].set_ylabel("exp_beta_ratio (log scale)")

axes[2].scatter(df["mean_tpm"], df["exp_beta_ratio"], s=2, alpha=0.2)
axes[2].set_xscale("log")
axes[2].set_yscale("log")
axes[2].set_xlabel("mean TPM across tissues (log scale)")
axes[2].set_ylabel("exp_beta_ratio (log scale)")

fig.suptitle(f"Gene specificity (tau) vs exp_beta exponential/beta likelihood ratio (n={len(df)} genes)")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)
