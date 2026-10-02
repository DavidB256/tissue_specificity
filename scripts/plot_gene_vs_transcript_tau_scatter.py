"""Scatter: one point per transcript of a multi-transcript gene (>=2 protein-coding
CHESS transcripts, per multitranscript_variance.csv), x=tau_unweighted of the
transcript's gene (compute_gene_tau_gtex_txs.py), y=tau_unweighted of the transcript
itself (compute_tau_unweighted.py). OLS fit y = b0 + b1*x (statsmodels), equation +
95% CIs for b0 and b1 reported in the legend."""
import sys
import numpy as np
import pandas as pd
import statsmodels.api as sm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IN_TX_TAU_CSV = sys.argv[1]
IN_META_CSV = sys.argv[2]
IN_GENE_TAU_CSV = sys.argv[3]
IN_MULTITX_CSV = sys.argv[4]
OUT_PNG = sys.argv[5]

tx_tau = pd.read_csv(IN_TX_TAU_CSV).rename(columns={"gene": "transcript_id", "tau_unweighted": "tau_tx"})
meta = pd.read_csv(IN_META_CSV)[["transcript_id", "gene_id"]]
gene_tau = pd.read_csv(IN_GENE_TAU_CSV)[["gene_id", "tau_unweighted"]].rename(columns={"tau_unweighted": "tau_gene"})
multitx_genes = pd.read_csv(IN_MULTITX_CSV)["gene_id"]

df = tx_tau.merge(meta, on="transcript_id", how="left").merge(gene_tau, on="gene_id", how="left")
df = df[df["gene_id"].isin(multitx_genes)]
df = df.dropna(subset=["tau_tx", "tau_gene"])

x = df["tau_gene"].to_numpy()
y = df["tau_tx"].to_numpy()
model = sm.OLS(y, sm.add_constant(x)).fit()
b0, b1 = model.params
(b0_lo, b0_hi), (b1_lo, b1_hi) = model.conf_int(alpha=0.05)

fig, ax = plt.subplots(figsize=(7, 7))
ax.scatter(x, y, s=8, alpha=0.1, color="steelblue")

x_line = np.array([x.min(), x.max()])
label = (
    f"y = {b0:.3f} + {b1:.3f} x\n"
    f"95% CI $\\beta_0$: [{b0_lo:.3f}, {b0_hi:.3f}]\n"
    f"95% CI $\\beta_1$: [{b1_lo:.3f}, {b1_hi:.3f}]"
)
ax.plot(x_line, b0 + b1 * x_line, color="red", linewidth=2, label=label)

ax.set_xlabel("gene tau_unweighted")
ax.set_ylabel("transcript tau_unweighted")
ax.set_title(f"gtex_txs: transcript tau vs. its gene's tau, multi-transcript genes only (n={len(df)} transcripts)")
ax.legend(loc="lower right", fontsize=9)
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)

print(f"n={len(df)} b0={b0:.4f} [{b0_lo:.4f},{b0_hi:.4f}] b1={b1:.4f} [{b1_lo:.4f},{b1_hi:.4f}] -> {OUT_PNG}",
      file=sys.stderr)
