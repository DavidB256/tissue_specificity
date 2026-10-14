"""Scatter: per gene, number of transcripts with isoform_fraction (tx_total /
gene_total TPM, summed across all tissues) >= 1% (x) vs. mean of
(tx_tau_weighted - gene_tau_weighted) across those transcripts (y). Points
colored by gene_tau_weighted. LOESS trend line (frac=0.2)."""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from statsmodels.nonparametric.smoothers_lowess import lowess

IN_GCT = sys.argv[1]
IN_TX_TAU_CSV = sys.argv[2]
IN_INFO_CONTENT_CSV = sys.argv[3]
OUT_PNG = sys.argv[4]

SIG_FRACTION = 0.01
LOESS_FRAC = 0.2

gct = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(gct.columns[2:])
tx_total = gct[tissues].sum(axis=1)
keep = tx_total > 0

tx = pd.DataFrame({
    "transcript_id": gct["Name"][keep],
    "gene_id": gct["Description"][keep],
    "tx_total": tx_total[keep],
})
tx["gene_total"] = tx.groupby("gene_id")["tx_total"].transform("sum")
tx["isoform_fraction"] = tx["tx_total"] / tx["gene_total"]

tx_tau = pd.read_csv(IN_TX_TAU_CSV).rename(columns={
    "gene": "transcript_id", "tau_free_weighted": "tx_tau_weighted",
})[["transcript_id", "tx_tau_weighted"]]
info = pd.read_csv(IN_INFO_CONTENT_CSV)[["gene_id", "gene_tau_weighted"]]

df = tx.merge(tx_tau, on="transcript_id", how="inner")
df = df.merge(info, on="gene_id", how="inner")
df = df.dropna(subset=["tx_tau_weighted", "gene_tau_weighted", "isoform_fraction"])
df["deviation"] = df["tx_tau_weighted"] - df["gene_tau_weighted"]

sig = df[df["isoform_fraction"] >= SIG_FRACTION]
gene_stats = sig.groupby("gene_id").agg(
    n_sig_tx=("transcript_id", "count"),
    mean_deviation=("deviation", "mean"),
    gene_tau_weighted=("gene_tau_weighted", "first"),
).reset_index()

fig, ax = plt.subplots(figsize=(10, 7))
sca = ax.scatter(gene_stats["n_sig_tx"], gene_stats["mean_deviation"], c=gene_stats["gene_tau_weighted"],
                  cmap="viridis", s=8, alpha=0.6)
fig.colorbar(sca, ax=ax, label="gene tau_weighted")

loess_fit = lowess(gene_stats["mean_deviation"], gene_stats["n_sig_tx"], frac=LOESS_FRAC)
ax.plot(loess_fit[:, 0], loess_fit[:, 1], color="red", linewidth=1.5)

ax.axhline(0, color="gray", linestyle="--", linewidth=1)
ax.set_xlabel(f"n tx with isoform_fraction >= {SIG_FRACTION:.0%}")
ax.set_ylabel("mean(tx tau_weighted - gene tau_weighted)")
ax.set_title(f"gtex_v11_tx_protein_coding: n significant tx vs. mean tau deviation (n_genes={len(gene_stats)})")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)

print(f"n_genes={len(gene_stats)} -> {OUT_PNG}", file=sys.stderr)
