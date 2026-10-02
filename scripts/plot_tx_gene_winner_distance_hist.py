"""Histogram of per-transcript winning_tissue_tx vs. winning_tissue_gene distance
(1 - Pearson r between the two winning tissues' profiles, from
compute_tx_gene_winner_distance.py). Explicitly reports the number of
transcripts won by the same tissue as their gene (dist=0)."""
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IN_CSV = sys.argv[1]
OUT_PNG = sys.argv[2]

df = pd.read_csv(IN_CSV)
n_total = len(df)
n_same = int((df["winning_tissue_tx"] == df["winning_tissue_gene"]).sum())

fig, ax = plt.subplots(figsize=(8, 6))
ax.hist(df["dist"], bins=50, color="steelblue", edgecolor="black", linewidth=0.3)
ax.set_xlabel("distance (1 - Pearson r) between tx's winning tissue and gene's winning tissue")
ax.set_ylabel("count of transcripts")
ax.set_title(
    f"tx vs. gene winning-tissue distance (n={n_total} transcripts)\n"
    f"same tissue as gene: {n_same} ({n_same / n_total * 100:.1f}%)"
)
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)

print(f"n={n_total} same_tissue_as_gene={n_same} ({n_same / n_total * 100:.1f}%) -> {OUT_PNG}", file=sys.stderr)
