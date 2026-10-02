"""Scatterplots, gene-level tau_unweighted (compute_gene_tau_gtex_txs.py) vs. (i)
tau_range and (ii) tau_var among that gene's own transcripts (compute_multitranscript_
variance.py -- multi-transcript genes only). Row 1 colors points by the gene's
median_log_tpm (median across the 9 tissues of log10(TPM+1), floor_tpm'd); row 2
repeats the same two (x, y) panels colored by n_transcripts (log-scaled -- heavily
right-skewed, median 4, max ~150) instead. A red LOESS curve (moderately wide
bandwidth, frac=0.4) is overlaid on every panel; points are drawn at alpha=0.1 so the
LOESS trend reads clearly through the dense overplotting."""
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from statsmodels.nonparametric.smoothers_lowess import lowess

LOESS_FRAC = 0.4

IN_GENE_TAU_CSV = sys.argv[1]
IN_MULTITX_CSV = sys.argv[2]
OUT_PNG = sys.argv[3]

gene_tau = pd.read_csv(IN_GENE_TAU_CSV)
multitx = pd.read_csv(IN_MULTITX_CSV)
df = gene_tau.merge(multitx, on="gene_id").dropna(subset=["tau_unweighted"])

fig, axes = plt.subplots(2, 2, figsize=(13, 11))
xy_panels = [("tau_range", "transcript tau_range"), ("tau_var", "transcript tau_var")]
color_rows = [
    ("median_log_tpm", "median_log_tpm (log10 TPM+1, median across tissues)", {}),
    ("n_transcripts", "n_transcripts (log scale)", {"norm": LogNorm()}),
]

for row, (color_col, color_label, extra_kwargs) in enumerate(color_rows):
    sc = None
    for col_idx, (col, ylabel) in enumerate(xy_panels):
        ax = axes[row, col_idx]
        sub = df.dropna(subset=[col])
        sc = ax.scatter(sub["tau_unweighted"], sub[col], c=sub[color_col],
                         cmap="viridis", s=10, alpha=0.1, **extra_kwargs)
        smoothed = lowess(sub[col], sub["tau_unweighted"], frac=LOESS_FRAC)
        ax.plot(smoothed[:, 0], smoothed[:, 1], color="red", linewidth=2)
        ax.set_xlabel("gene tau_unweighted")
        ax.set_ylabel(ylabel)
        ax.set_title(f"{ylabel} vs gene tau (n={len(sub)} genes)")
    fig.colorbar(sc, ax=axes[row, :], label=color_label)

fig.suptitle("gtex_txs: gene tau vs. spread among its own transcripts' tau")
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)

print(f"genes={len(df)} -> {OUT_PNG}", file=sys.stderr)
