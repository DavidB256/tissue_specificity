"""Scatter comparing a per-gene paired variance (default: delete-one-cluster jackknife
variance, see jackknife_cluster_variance.py; also reused for Monte Carlo
noise-perturbation variance, see compute_noise_variance.py) between tau_unweighted and
a weighted variant: one point per gene, x=unweighted variance, y=weighted variance,
colored by that gene's tau_weighted point estimate (the paired mean column -- same
name as WEIGHTED_COL with "_var" swapped for "_mean", already present in both summary
CSV families). Wilcoxon signed-rank test (paired, non-parametric -- variances are
strongly right-skewed) on the per-gene variance difference -- this is the
null-hypothesis test that both estimators undergo equal variation under the same
perturbation."""
import sys
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plot_utils import display_label, annotate_diag_stats

IN_SUMMARY = sys.argv[1]
OUT_PNG = sys.argv[2]
LABEL = sys.argv[3] if len(sys.argv) > 3 else "reference"
UNWEIGHTED_COL = sys.argv[4] if len(sys.argv) > 4 else "tau_unweighted_jk_var"
WEIGHTED_COL = sys.argv[5] if len(sys.argv) > 5 else "tau_free_cluster_weighted_jk_var"
VAR_KIND = sys.argv[6] if len(sys.argv) > 6 else "jackknife"
WEIGHTED_NAME = sys.argv[7] if len(sys.argv) > 7 else "free_cluster_weighted"
COLOR_COL = WEIGHTED_COL.replace("_var", "_mean")

df = pd.read_csv(IN_SUMMARY)[["gene", UNWEIGHTED_COL, WEIGHTED_COL, COLOR_COL]].dropna()
u = df[UNWEIGHTED_COL].values
w = df[WEIGHTED_COL].values
c = df[COLOR_COL].values

stat, pval = wilcoxon(u, w)
median_u, median_w = float(np.median(u)), float(np.median(w))
mean_u, mean_w = float(u.mean()), float(w.mean())

# floor well above the near-zero pile from genes whose tau is algebraically invariant to
# reweighting (e.g. single-tissue-specific genes, tau=1 regardless of w) -- those are a
# real but separate phenomenon from the continuous variance distribution being compared,
# and left unfloored they'd stretch the log axis into an uninformative multi-decade spike
eps = 1e-6

log_u = np.log10(np.clip(u, eps, None))
log_w = np.log10(np.clip(w, eps, None))

fig, ax = plt.subplots(figsize=(8.5, 7))

sc = ax.scatter(log_u, log_w, c=c, s=6, alpha=0.6, cmap="viridis", rasterized=True)
cb = fig.colorbar(sc, ax=ax)
cb.set_label(f"tau_{WEIGHTED_NAME} (point estimate)")

lo = min(log_u.min(), log_w.min())
hi = max(log_u.max(), log_w.max())
ax.plot([lo, hi], [lo, hi], "--", color="gray", lw=1, zorder=1)
ax.set_xlim(lo, hi)
ax.set_ylim(lo, hi)
ax.set_xlabel(f"log10(tau_unweighted {VAR_KIND} variance) [floored at {eps:g}]")
ax.set_ylabel(f"log10(tau_{WEIGHTED_NAME} {VAR_KIND} variance) [floored at {eps:g}]")
r, rho, frac_below, frac_above = annotate_diag_stats(ax, log_u, log_w)
ax.set_title(f"{display_label(LABEL)}: {VAR_KIND} variance, unweighted vs {WEIGHTED_NAME} "
             f"(n={len(df)} genes)\nWilcoxon signed-rank p={pval:.3g}, "
             f"median {median_u:.2e} -> {median_w:.2e}", fontsize=10)
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)

print(f"genes={len(df)} median_var_unweighted={median_u:.6g} median_var_{WEIGHTED_NAME}={median_w:.6g} "
      f"mean_var_unweighted={mean_u:.6g} mean_var_{WEIGHTED_NAME}={mean_w:.6g} "
      f"wilcoxon_stat={stat:.6g} pval={pval:.6g} frac_below={frac_below:.3f} frac_above={frac_above:.3f}",
      file=sys.stderr)
print(f"wrote {OUT_PNG}", file=sys.stderr)
