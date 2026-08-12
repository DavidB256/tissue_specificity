"""Same as plot_tau_bootstrap_variance.py, restricted to high-specificity genes:
tau_unweighted >= 0.8 OR tau_weighted_gls_ridge >= 0.8 (canonical values from
weighted_tau.csv, not the bootstrap mean, so the same gene set is used across all
three resampling panels)."""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plot_utils import bootstrap_scatter_panel

IN_TAU = sys.argv[1]
IN_UNWEIGHTED_TISSUE = sys.argv[2]
IN_RIDGE_TISSUE = sys.argv[3]
IN_UNWEIGHTED_CLUSTER10 = sys.argv[4]
IN_RIDGE_CLUSTER10 = sys.argv[5]
IN_UNWEIGHTED_CLUSTER20 = sys.argv[6]
IN_RIDGE_CLUSTER20 = sys.argv[7]
OUT_PNG = sys.argv[8]

VAR_FLOOR = 1e-5
TAU_THRESHOLD = 0.8

tau = pd.read_csv(IN_TAU)
high_spec_genes = set(tau.loc[
    (tau["tau_unweighted"] >= TAU_THRESHOLD) | (tau["tau_weighted_gls_ridge"] >= TAU_THRESHOLD),
    "gene"
])
print(f"high-specificity genes (tau_unweighted>={TAU_THRESHOLD} or "
      f"tau_weighted_gls_ridge>={TAU_THRESHOLD}): {len(high_spec_genes)}/{len(tau)}", file=sys.stderr)


def load_pair(unweighted_csv, ridge_csv, mode):
    ucol, rcol = f"tau_unweighted_{mode}_var", f"tau_ridge_{mode}_var"
    mcol = f"tau_unweighted_{mode}_mean"
    u = pd.read_csv(unweighted_csv)[["gene", ucol, mcol]]
    w = pd.read_csv(ridge_csv)[["gene", rcol]]
    df = u.merge(w, on="gene", how="inner").replace([np.inf, -np.inf], np.nan).dropna()
    df = df[df["gene"].isin(high_spec_genes)]
    return df, ucol, rcol, mcol


fig, axes = plt.subplots(1, 3, figsize=(19, 5.5))

for ax, unweighted_csv, ridge_csv, mode, title in [
    (axes[0], IN_UNWEIGHTED_TISSUE, IN_RIDGE_TISSUE, "tissue", "tissue resampling"),
    (axes[1], IN_UNWEIGHTED_CLUSTER10, IN_RIDGE_CLUSTER10, "cluster10", "correlation-block resampling (t=10)"),
    (axes[2], IN_UNWEIGHTED_CLUSTER20, IN_RIDGE_CLUSTER20, "cluster20", "correlation-block resampling (t=20)"),
]:
    df, ucol, rcol, mcol = load_pair(unweighted_csv, ridge_csv, mode)
    bootstrap_scatter_panel(
        fig, ax, df[ucol], df[rcol], df[mcol],
        title=f"{title} (n={len(df)} genes)",
        xlabel="Var(tau_unweighted) across resamples",
        ylabel="Var(tau_weighted_gls_ridge) across resamples",
        floor=VAR_FLOOR,
    )
    print(f"{mode}: n={len(df)}", file=sys.stderr)

fig.suptitle(f"Bootstrap variance of tau, high-specificity genes "
             f"(tau_unweighted>={TAU_THRESHOLD} or tau_weighted_gls_ridge>={TAU_THRESHOLD})")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)
