"""Scatter: one point per transcript, x=n_transcripts (protein-coding transcript
count of the tx's gene), y=tx-gene winner distance (1 - Pearson r between the
tx's winning tissue and its gene's winning tissue,
compute_tx_gene_winner_distance.py). Points colored by gene tau_unweighted.
Red LOESS curve (moderate bandwidth, frac=LOESS_FRAC) overlaid."""
import sys
import numpy as np
import pandas as pd
from statsmodels.nonparametric.smoothers_lowess import lowess
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

LOESS_FRAC = 0.4

IN_DIST_CSV = sys.argv[1]
IN_GENE_TAU_CSV = sys.argv[2]
IN_GCT = sys.argv[3]
OUT_PNG = sys.argv[4]
TXS_SOURCE = sys.argv[5]

dist = pd.read_csv(IN_DIST_CSV)[["transcript_id", "gene_id", "dist"]]
gene_tau = pd.read_csv(IN_GENE_TAU_CSV)[["gene_id", "tau_unweighted"]].rename(columns={"tau_unweighted": "tau_gene"})

gct = pd.read_csv(IN_GCT, sep="\t", skiprows=2, usecols=["Description"])
n_transcripts = gct.groupby("Description").size().rename("n_transcripts")

df = dist.merge(gene_tau, on="gene_id", how="left").merge(n_transcripts, left_on="gene_id", right_index=True, how="left")
df = df.dropna(subset=["dist", "n_transcripts", "tau_gene"])

fig, ax = plt.subplots(figsize=(8, 7))
sc = ax.scatter(df["n_transcripts"], df["dist"], c=df["tau_gene"], cmap="viridis", s=8, alpha=0.1)
smoothed = lowess(df["dist"], df["n_transcripts"], frac=LOESS_FRAC)
ax.plot(smoothed[:, 0], smoothed[:, 1], color="red", linewidth=2)
ax.set_xlabel("n_transcripts (of tx's gene)")
ax.set_ylabel("tx-gene winner distance (1 - Pearson r)")
ax.set_title(f"{TXS_SOURCE}: tx-gene winner distance vs. n_transcripts (n={len(df)} transcripts)")
fig.colorbar(sc, ax=ax, label="gene tau_unweighted")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)

print(f"n={len(df)} -> {OUT_PNG}", file=sys.stderr)
