"""Scatter: one point per multi-transcript gene (n_tx>=2 -- single-tx genes are
trivially weighted_kl=0, uninformative here), x=weighted_kl (TPM-share-weighted
average of KL(tx tissue pmf || gene tissue pmf), bits), y=gene_tau_weighted
(gene's own tau_free_weighted, hierclust0.1/max-agg) -- both from
compute_tx_info_content.py's tx_info_content.csv. OLS fit (equation + 95% CIs
for b0/b1 in the legend) and a red LOESS curve overlaid."""
import sys
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.nonparametric.smoothers_lowess import lowess
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

LOESS_FRAC = 0.4

IN_CSV = sys.argv[1]
OUT_PNG = sys.argv[2]

df = pd.read_csv(IN_CSV).dropna(subset=["weighted_kl", "gene_tau_weighted"])
df = df[df["n_tx"] >= 2]
x = df["weighted_kl"].to_numpy()
y = df["gene_tau_weighted"].to_numpy()

model = sm.OLS(y, sm.add_constant(x)).fit()
b0, b1 = model.params
(b0_lo, b0_hi), (b1_lo, b1_hi) = model.conf_int(alpha=0.05)

fig, ax = plt.subplots(figsize=(8, 7))
ax.scatter(x, y, s=8, alpha=0.1, color="steelblue")

x_line = np.array([x.min(), x.max()])
label = (
    f"y = {b0:.3f} + {b1:.3f} x\n"
    f"95% CI $\\beta_0$: [{b0_lo:.3f}, {b0_hi:.3f}]\n"
    f"95% CI $\\beta_1$: [{b1_lo:.3f}, {b1_hi:.3f}]"
)
ax.plot(x_line, b0 + b1 * x_line, color="red", linewidth=2, label=label)

smoothed = lowess(y, x, frac=LOESS_FRAC)
ax.plot(smoothed[:, 0], smoothed[:, 1], color="darkorange", linewidth=2, label="LOESS")

ax.set_xlabel("weighted KL sum (bits)")
ax.set_ylabel("gene tau_weighted (free_cluster_weights, hierclust0.1, max-agg)")
ax.set_title(f"gtex_v11_tx_protein_coding: weighted KL vs. gene tau_weighted, multi-tx genes only (n={len(df)} genes)")
ax.legend(loc="lower right", fontsize=9)
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)

print(f"n={len(df)} b0={b0:.4f} [{b0_lo:.4f},{b0_hi:.4f}] b1={b1:.4f} [{b1_lo:.4f},{b1_hi:.4f}] -> {OUT_PNG}",
      file=sys.stderr)
