"""Visualize per-gene bootstrap variance of tau_unweighted vs tau_weighted_gls_ridge,
under (i) tissue resampling, (ii) t=10 correlation-block resampling, and (iii) t=20
correlation-block resampling. Log-log scatter since variances span orders of magnitude
across genes; displayed variances are left-censored at 1e-5."""
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

VAR_FLOOR = 1e-5  # left-censor displayed variances at this floor


def load_pair(unweighted_csv, ridge_csv, mode):
    ucol, rcol = f"tau_unweighted_{mode}_var", f"tau_ridge_{mode}_var"
    mcol = f"tau_unweighted_{mode}_mean"
    u = pd.read_csv(unweighted_csv)[["gene", ucol, mcol]]
    w = pd.read_csv(ridge_csv)[["gene", rcol]]
    df = u.merge(w, on="gene", how="inner").replace([np.inf, -np.inf], np.nan).dropna()
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

fig.suptitle("Bootstrap variance of tau: unweighted vs ridge-weighted")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)
