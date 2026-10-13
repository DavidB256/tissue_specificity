"""Histogram of arctanh(tx_tau_weighted) - arctanh(gene_tau_weighted) across
all quantified transcripts (tx_total>0) of multi-tx genes (n_tx>=2). Each tx
weighted by 1/n_tx of its gene, so every gene contributes equal total mass
(a 5-tx gene's transcripts each weigh 1/5) regardless of isoform count.
tau scaled by 0.9999 before arctanh (arctanh(1) undefined)."""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IN_GCT = sys.argv[1]
IN_TX_TAU_CSV = sys.argv[2]
IN_INFO_CONTENT_CSV = sys.argv[3]
OUT_PNG = sys.argv[4]

SCALE = 0.9999


def arctanh_tau(tau):
    return np.arctanh(SCALE * tau)


gct = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(gct.columns[2:])
tx_total = gct[tissues].sum(axis=1)
keep = tx_total > 0

tx = pd.DataFrame({"transcript_id": gct["Name"][keep], "gene_id": gct["Description"][keep]})

tx_tau = pd.read_csv(IN_TX_TAU_CSV).rename(columns={
    "gene": "transcript_id", "tau_free_weighted": "tx_tau_weighted",
})[["transcript_id", "tx_tau_weighted"]]
info = pd.read_csv(IN_INFO_CONTENT_CSV)[["gene_id", "gene_tau_weighted"]]

df = tx.merge(tx_tau, on="transcript_id", how="inner")
df = df.merge(info, on="gene_id", how="inner")
df = df.dropna(subset=["tx_tau_weighted", "gene_tau_weighted"])

n_tx_per_gene = df.groupby("gene_id")["transcript_id"].transform("count")
df = df[n_tx_per_gene >= 2].copy()
n_tx_per_gene = n_tx_per_gene[df.index]

df["deviation_arctanh"] = arctanh_tau(df["tx_tau_weighted"]) - arctanh_tau(df["gene_tau_weighted"])
df["weight"] = 1.0 / n_tx_per_gene

fig, ax = plt.subplots(figsize=(10, 7))
ax.hist(df["deviation_arctanh"], weights=df["weight"], bins=100, color="steelblue", alpha=0.7)
ax.axvline(0, color="gray", linestyle="--", linewidth=1)
ax.set_xlabel("arctanh(tx tau_weighted) - arctanh(gene tau_weighted)")
ax.set_ylabel("gene-weighted count (tx weighted by 1/n_tx of its gene)")
ax.set_title(
    f"gtex_v11_tx_protein_coding: arctanh tau deviation, gene-weighted\n"
    f"n_tx={len(df)}, n_genes={df['gene_id'].nunique()}, total weight={df['weight'].sum():.0f}"
)
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)

print(f"n_tx={len(df)} n_genes={df['gene_id'].nunique()} total_weight={df['weight'].sum():.4f} -> {OUT_PNG}",
      file=sys.stderr)
