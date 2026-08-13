"""Unweighted tau under the full 68-tissue panel vs 'neural-downsampled' tau (mean tau
across the 13 subsets of all-non-brain-tissues-plus-one-brain-tissue; see
compute_tau_neural_downsample.py). Isolates how much brain's 13x subregion redundancy
inflates tau_unweighted, with no GLS/ridge weighting involved. Genes restricted to
exactly 1 non-brain tissue at TPM>=1 are excluded: dropping brain to 1 tissue trivially
dominates their tau in the downsampled subset."""
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plot_utils import annotate_diag_stats

IN_TAU = sys.argv[1]
IN_DOWNSAMPLE = sys.argv[2]
OUT_PNG = sys.argv[3]

tau = pd.read_csv(IN_TAU)[["gene", "tau_unweighted"]]
ds = pd.read_csv(IN_DOWNSAMPLE)[["gene", "n_nonbrain_tpm1", "tau_neural_downsampled_mean"]]
df = tau.merge(ds, on="gene", how="inner")

n_before = len(df)
df = df[df["n_nonbrain_tpm1"] != 1].dropna(subset=["tau_unweighted", "tau_neural_downsampled_mean"])
print(f"excluded {n_before - len(df)}/{n_before} genes (restricted to 1 non-brain tissue, "
      f"or NaN tau); n={len(df)} remain", file=sys.stderr)

df["delta_tau"] = df["tau_neural_downsampled_mean"] - df["tau_unweighted"]

fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

axes[0].hist(df["delta_tau"], bins=100, color="tab:blue")
axes[0].axvline(0, ls="--", color="gray", lw=1)
axes[0].set_xlabel("delta_tau = tau_neural_downsampled_mean - tau_unweighted")
axes[0].set_ylabel("count")
axes[0].set_title(f"delta_tau distribution (n={len(df)} genes)", fontsize=10)

axes[1].scatter(df["tau_unweighted"], df["tau_neural_downsampled_mean"], s=4, alpha=0.15, color="tab:blue")
axes[1].plot([0, 1], [0, 1], "--", color="gray", lw=1)
axes[1].set_xlabel("tau_unweighted (full 68-tissue panel)")
axes[1].set_ylabel("tau_neural_downsampled_mean (13x 1-brain-tissue subset mean)")
axes[1].set_title(f"tau_unweighted vs neural-downsampled (n={len(df)} genes)", fontsize=10)
annotate_diag_stats(axes[1], df["tau_unweighted"], df["tau_neural_downsampled_mean"])

fig.suptitle("Effect of downsampling brain to 1 (of 13) tissues on unweighted tau")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)
print(f"delta_tau: mean={df['delta_tau'].mean():.4f} median={df['delta_tau'].median():.4f}", file=sys.stderr)
