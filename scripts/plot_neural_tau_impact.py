"""Impact of GLS/ridge tau weighting specifically on 'neural' genes (max expression
in a Brain_* tissue). The 13 GTEx brain subregions are highly correlated, so this is
the case study where the weighting should matter most: a gene expressed similarly
across several brain regions gets an inflated unweighted tau (each correlated region
counted as independent evidence of specificity), which GLS/ridge weighting should
shrink."""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IN_TAU = sys.argv[1]
IN_NEURAL = sys.argv[2]
OUT_CSV = sys.argv[3]
OUT_PNG = sys.argv[4]

tau = pd.read_csv(IN_TAU)
neural = pd.read_csv(IN_NEURAL)[["gene", "max_tissue", "is_neural"]]
df = tau.merge(neural, on="gene", how="inner").dropna(
    subset=["tau_unweighted", "tau_weighted_gls", "tau_weighted_gls_ridge"])

df["delta_gls"] = df["tau_weighted_gls"] - df["tau_unweighted"]
df["delta_gls_ridge"] = df["tau_weighted_gls_ridge"] - df["tau_unweighted"]

neural_df = df[df.is_neural].copy()
neural_df[["gene", "max_tissue", "tau_unweighted", "tau_weighted_gls", "tau_weighted_gls_ridge",
           "delta_gls", "delta_gls_ridge"]].to_csv(OUT_CSV, index=False)

fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))

for ax, ycol, label, robust in [(axes[0], "tau_weighted_gls", "GLS weighted (raw)", True),
                                 (axes[1], "tau_weighted_gls_ridge", "GLS ridge weighted", False)]:
    ax.scatter(neural_df["tau_unweighted"], neural_df[ycol], s=4, alpha=0.2, color="tab:blue")
    ax.plot([0, 1], [0, 1], "--", color="gray", lw=1)
    ax.set_xlim(0, 1)
    ax.set_xlabel("tau_unweighted")
    ax.set_ylabel(ycol)
    if robust:
        lo, hi = np.nanpercentile(neural_df[ycol], [1, 99])
        pad = 0.05 * (hi - lo)
        ax.set_ylim(lo - pad, hi + pad)
        ax.set_title(f"neural genes (n={len(neural_df)}): unweighted vs {label}\n(y clipped to 1st-99th pctile)", fontsize=9)
    else:
        ax.set_ylim(0, 1)
        ax.set_title(f"neural genes (n={len(neural_df)}): unweighted vs {label}", fontsize=10)

ax = axes[2]
bins = np.linspace(-1, 1, 61)
ax.hist(df["delta_gls_ridge"], bins=bins, density=True, alpha=0.5, label=f"all genes (n={len(df)})", color="gray")
ax.hist(neural_df["delta_gls_ridge"], bins=bins, density=True, alpha=0.5, label=f"neural genes (n={len(neural_df)})", color="tab:red")
ax.axvline(0, color="black", lw=1)
ax.set_xlabel("delta = tau_weighted_gls_ridge - tau_unweighted")
ax.set_ylabel("density")
ax.legend(fontsize=8)
ax.set_title("weighting impact: neural vs all genes", fontsize=10)

fig.suptitle("Impact of GLS/ridge tau weighting on neural genes (max expression in a Brain_* tissue)")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)

print(f"neural genes: {len(neural_df)}; median delta_gls_ridge neural={neural_df['delta_gls_ridge'].median():.4f} "
      f"vs all={df['delta_gls_ridge'].median():.4f}", file=sys.stderr)
