"""Visualize per-gene bootstrap coefficient of variation (CV = variance / mean tau,
as specified) of tau_unweighted vs tau_weighted_gls_ridge, under (i) tissue resampling,
(ii) t=10 correlation-block resampling, and (iii) t=20 correlation-block resampling.
Parallel figure to plot_tau_bootstrap_variance.py, using the same variance CSVs.
Genes with near-zero mean tau are dropped (CV undefined/explosive there); displayed CV
is left-censored at 1e-5, matching the variance figure's convention."""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plot_utils import bootstrap_scatter_panel

IN_UNWEIGHTED_TISSUE = sys.argv[1]
IN_RIDGE_TISSUE = sys.argv[2]
IN_UNWEIGHTED_CLUSTER10 = sys.argv[3]
IN_RIDGE_CLUSTER10 = sys.argv[4]
IN_UNWEIGHTED_CLUSTER20 = sys.argv[5]
IN_RIDGE_CLUSTER20 = sys.argv[6]
OUT_PNG = sys.argv[7]

CV_FLOOR = 1e-5
MEAN_MIN = 1e-3  # drop genes whose mean tau is too close to 0 for CV to be meaningful


def load_pair(unweighted_csv, ridge_csv, mode):
    uvar, umean = f"tau_unweighted_{mode}_var", f"tau_unweighted_{mode}_mean"
    rvar, rmean = f"tau_ridge_{mode}_var", f"tau_ridge_{mode}_mean"
    u = pd.read_csv(unweighted_csv)[["gene", uvar, umean]]
    w = pd.read_csv(ridge_csv)[["gene", rvar, rmean]]
    df = u.merge(w, on="gene", how="inner").replace([np.inf, -np.inf], np.nan).dropna()
    df = df[(df[umean].abs() > MEAN_MIN) & (df[rmean].abs() > MEAN_MIN)]
    df["cv_unweighted"] = df[uvar] / df[umean]
    df["cv_ridge"] = df[rvar] / df[rmean]
    return df


fig, axes = plt.subplots(1, 3, figsize=(19, 5.5))

for ax, unweighted_csv, ridge_csv, mode, title in [
    (axes[0], IN_UNWEIGHTED_TISSUE, IN_RIDGE_TISSUE, "tissue", "tissue resampling"),
    (axes[1], IN_UNWEIGHTED_CLUSTER10, IN_RIDGE_CLUSTER10, "cluster10", "correlation-block resampling (t=10)"),
    (axes[2], IN_UNWEIGHTED_CLUSTER20, IN_RIDGE_CLUSTER20, "cluster20", "correlation-block resampling (t=20)"),
]:
    df = load_pair(unweighted_csv, ridge_csv, mode)
    bootstrap_scatter_panel(
        fig, ax, df["cv_unweighted"], df["cv_ridge"], df[f"tau_unweighted_{mode}_mean"],
        title=f"{title} (n={len(df)} genes)",
        xlabel="CV(tau_unweighted) = Var/mean across resamples",
        ylabel="CV(tau_weighted_gls_ridge) = Var/mean across resamples",
        floor=CV_FLOOR,
    )
    print(f"{mode}: n={len(df)}", file=sys.stderr)

fig.suptitle("Bootstrap coefficient of variation of tau: unweighted vs ridge-weighted")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)
