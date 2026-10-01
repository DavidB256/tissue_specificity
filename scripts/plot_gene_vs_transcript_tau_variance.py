"""Two scatterplots, gene-level tau_unweighted (compute_gene_tau_gtex_txs.py) vs. (i)
tau_range and (ii) tau_var among that gene's own transcripts (compute_multitranscript_
variance.py -- multi-transcript genes only), points colored by the gene's
median_log_tpm (median across the 9 tissues of log10(TPM+1), floor_tpm'd)."""
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IN_GENE_TAU_CSV = sys.argv[1]
IN_MULTITX_CSV = sys.argv[2]
OUT_PNG = sys.argv[3]

gene_tau = pd.read_csv(IN_GENE_TAU_CSV)
multitx = pd.read_csv(IN_MULTITX_CSV)
df = gene_tau.merge(multitx, on="gene_id").dropna(subset=["tau_unweighted"])

fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
panels = [("tau_range", "transcript tau_range"), ("tau_var", "transcript tau_var")]
sc = None
for ax, (col, ylabel) in zip(axes, panels):
    sub = df.dropna(subset=[col])
    sc = ax.scatter(sub["tau_unweighted"], sub[col], c=sub["median_log_tpm"],
                     cmap="viridis", s=10, alpha=0.6)
    ax.set_xlabel("gene tau_unweighted")
    ax.set_ylabel(ylabel)
    ax.set_title(f"{ylabel} vs gene tau (n={len(sub)} genes)")

fig.colorbar(sc, ax=axes, label="median_log_tpm (log10 TPM+1, median across tissues)")
fig.suptitle("gtex_txs: gene tau vs. spread among its own transcripts' tau")
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)

print(f"genes={len(df)} -> {OUT_PNG}", file=sys.stderr)
