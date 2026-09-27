"""Compares the distributions of a per-gene paired variance (default: delete-one-
cluster jackknife variance, see jackknife_cluster_variance.py; also reused for Monte
Carlo noise-perturbation variance, see compute_noise_variance.py) between
tau_unweighted and a weighted variant, paired per gene. Wilcoxon signed-rank test
(paired, non-parametric -- variances are strongly right-skewed) on the per-gene
variance difference -- this is the null-hypothesis test that both estimators undergo
equal variation under the same perturbation."""
import sys
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from plot_utils import display_label

IN_SUMMARY = sys.argv[1]
OUT_PNG = sys.argv[2]
LABEL = sys.argv[3] if len(sys.argv) > 3 else "reference"
UNWEIGHTED_COL = sys.argv[4] if len(sys.argv) > 4 else "tau_unweighted_jk_var"
WEIGHTED_COL = sys.argv[5] if len(sys.argv) > 5 else "tau_free_cluster_weighted_jk_var"
VAR_KIND = sys.argv[6] if len(sys.argv) > 6 else "jackknife"
WEIGHTED_NAME = sys.argv[7] if len(sys.argv) > 7 else "free_cluster_weighted"

df = pd.read_csv(IN_SUMMARY)[["gene", UNWEIGHTED_COL, WEIGHTED_COL]].dropna()
u = df[UNWEIGHTED_COL].values
w = df[WEIGHTED_COL].values

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
n_up, n_down = int((w > u).sum()), int((w < u).sum())

fig, ax = plt.subplots(figsize=(9, 6.5))

segments = np.stack([np.column_stack([np.full(len(df), 1.0), log_u]),
                      np.column_stack([np.full(len(df), 2.0), log_w])], axis=1)
colors = np.where(w > u, "tab:red", np.where(w < u, "tab:blue", "tab:gray"))
lc = LineCollection(segments, colors=colors, linewidths=0.3, alpha=0.03, zorder=1)
ax.add_collection(lc)

parts = ax.violinplot([log_u, log_w], showmedians=True)
for body in parts["bodies"]:
    body.set_alpha(0.3)
    body.set_zorder(2)
for key in ("cmins", "cmaxes", "cbars", "cmedians"):
    parts[key].set_zorder(2)

ax.set_xticks([1, 2])
ax.set_xticklabels(["tau_unweighted", f"tau_{WEIGHTED_NAME}"])
ax.set_xlim(0.5, 2.5)
ax.set_ylim(-7, max(np.log10(max(u.max(), w.max())), 0) + 0.5)
ax.set_ylabel(f"log10({VAR_KIND} variance) [floored at {eps:g}]")
ax.set_title(f"{display_label(LABEL)}: {VAR_KIND} variance, unweighted vs {WEIGHTED_NAME} "
             f"(n={len(df)} genes)\nWilcoxon signed-rank p={pval:.3g}, "
             f"median {median_u:.2e} -> {median_w:.2e}\n"
             f"per-gene: {n_up} increase (red), {n_down} decrease (blue)", fontsize=10)
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)

print(f"genes={len(df)} median_var_unweighted={median_u:.6g} median_var_{WEIGHTED_NAME}={median_w:.6g} "
      f"mean_var_unweighted={mean_u:.6g} mean_var_{WEIGHTED_NAME}={mean_w:.6g} "
      f"wilcoxon_stat={stat:.6g} pval={pval:.6g}", file=sys.stderr)
print(f"wrote {OUT_PNG}", file=sys.stderr)
