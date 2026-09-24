"""Compares the distributions of delete-one-cluster jackknife variance (see
jackknife_cluster_variance.py) between tau_unweighted and tau_free_cluster_weighted,
paired per gene. Wilcoxon signed-rank test (paired, non-parametric -- variances are
strongly right-skewed) on the per-gene variance difference."""
import sys
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plot_utils import display_label

IN_JACKKNIFE_SUMMARY = sys.argv[1]
OUT_PNG = sys.argv[2]
LABEL = sys.argv[3] if len(sys.argv) > 3 else "reference"

df = pd.read_csv(IN_JACKKNIFE_SUMMARY)[["gene", "tau_unweighted_jk_var", "tau_free_cluster_weighted_jk_var"]].dropna()
u = df["tau_unweighted_jk_var"].values
w = df["tau_free_cluster_weighted_jk_var"].values

stat, pval = wilcoxon(u, w)
median_u, median_w = float(np.median(u)), float(np.median(w))
mean_u, mean_w = float(u.mean()), float(w.mean())

# floor well above the near-zero pile from genes whose tau is algebraically invariant to
# reweighting (e.g. single-tissue-specific genes, tau=1 regardless of w) -- those are a
# real but separate phenomenon from the continuous variance distribution being compared,
# and left unfloored they'd stretch the log axis into an uninformative multi-decade spike
eps = 1e-6

fig, ax = plt.subplots(figsize=(9, 6.5))
ax.violinplot([np.log10(np.clip(u, eps, None)), np.log10(np.clip(w, eps, None))], showmedians=True)
ax.set_xticks([1, 2])
ax.set_xticklabels(["tau_unweighted", "tau_free_cluster_weighted"])
ax.set_ylim(-7, max(np.log10(max(u.max(), w.max())), 0) + 0.5)
ax.set_ylabel(f"log10(jackknife variance) [floored at {eps:g}]")
ax.set_title(f"{display_label(LABEL)}: jackknife variance, unweighted vs free_cluster_weighted "
             f"(n={len(df)} genes)\nWilcoxon signed-rank p={pval:.3g}, "
             f"median {median_u:.2e} -> {median_w:.2e}", fontsize=10)
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)

print(f"genes={len(df)} median_var_unweighted={median_u:.6g} median_var_free_cluster_weighted={median_w:.6g} "
      f"mean_var_unweighted={mean_u:.6g} mean_var_free_cluster_weighted={mean_w:.6g} "
      f"wilcoxon_stat={stat:.6g} pval={pval:.6g}", file=sys.stderr)
print(f"wrote {OUT_PNG}", file=sys.stderr)
