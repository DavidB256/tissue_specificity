"""tau vs Var(tau) scatter grid: 2 rows (tau_unweighted, tau_weighted_gls_ridge) x 3
cols (tissue, cluster10, cluster20 resampling). Shows whether bootstrap variance of a
gene's tau depends on the gene's tau value itself, for each tau variant/resampling
scheme. Displayed variance is left-censored at 1e-5, matching the other bootstrap
figures."""
import sys
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IN_UNWEIGHTED_TISSUE = sys.argv[1]
IN_RIDGE_TISSUE = sys.argv[2]
IN_UNWEIGHTED_CLUSTER10 = sys.argv[3]
IN_RIDGE_CLUSTER10 = sys.argv[4]
IN_UNWEIGHTED_CLUSTER20 = sys.argv[5]
IN_RIDGE_CLUSTER20 = sys.argv[6]
OUT_PNG = sys.argv[7]

VAR_FLOOR = 1e-5

MODES = [
    ("tissue", IN_UNWEIGHTED_TISSUE, IN_RIDGE_TISSUE, "tissue resampling"),
    ("cluster10", IN_UNWEIGHTED_CLUSTER10, IN_RIDGE_CLUSTER10, "correlation-block resampling (t=10)"),
    ("cluster20", IN_UNWEIGHTED_CLUSTER20, IN_RIDGE_CLUSTER20, "correlation-block resampling (t=20)"),
]

TAU_ROWS = [("unweighted", "tau_unweighted"), ("ridge", "tau_weighted_gls_ridge")]

fig, axes = plt.subplots(2, 3, figsize=(19, 10.5))

for row, (tau_key, tau_label) in enumerate(TAU_ROWS):
    for col, (mode, unweighted_csv, ridge_csv, title) in enumerate(MODES):
        csv = unweighted_csv if tau_key == "unweighted" else ridge_csv
        mcol, vcol = f"tau_{tau_key}_{mode}_mean", f"tau_{tau_key}_{mode}_var"
        df = pd.read_csv(csv)[["gene", mcol, vcol]].replace([np.inf, -np.inf], np.nan).dropna()
        x = df[mcol]
        y = np.maximum(df[vcol], VAR_FLOOR)

        ax = axes[row, col]
        ax.scatter(x, y, s=4, alpha=0.1, color="tab:blue")
        ax.set_yscale("log")
        ax.set_xlim(0, 1)
        ax.set_xlabel(tau_label)
        ax.set_ylabel(f"Var({tau_label})")
        r, _ = pearsonr(x, y)
        rho, _ = spearmanr(x, y)
        ax.text(0.03, 0.97, f"pearson r={r:.3f}\nspearman ρ={rho:.3f}",
                transform=ax.transAxes, va="top", ha="left", fontsize=9,
                bbox=dict(facecolor="white", alpha=0.7, edgecolor="none"))
        ax.set_title(f"{title} (n={len(df)} genes)", fontsize=10)
        print(f"{tau_key}/{mode}: n={len(df)} pearson={r:.4f} spearman={rho:.4f}", file=sys.stderr)

fig.suptitle("Gene tau value vs bootstrap variance of tau")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)
