"""Histogram of arctanh(tx_tau_weighted) - arctanh(gene_tau_weighted) across
multi-tx genes (n_tx>=2, counted within each panel's own filter). Each tx
weighted by 1/n_tx of its gene (within that panel), so every gene
contributes equal total mass (a 5-tx gene's transcripts each weigh 1/5)
regardless of isoform count. tau scaled by 0.9999 before arctanh
(arctanh(1) undefined). Panel 1: all quantified transcripts (tx_total>0).
Panel 2: only tx with isoform_fraction (tx_total/gene_total TPM, summed
across all tissues) >= 1%."""
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

base = tx.merge(tx_tau, on="transcript_id", how="inner")
base = base.merge(info, on="gene_id", how="inner")
base = base.dropna(subset=["tx_tau_weighted", "gene_tau_weighted", "isoform_fraction"])
base["deviation_arctanh"] = arctanh_tau(base["tx_tau_weighted"]) - arctanh_tau(base["gene_tau_weighted"])

PANELS = [("all quantified tx", 0.0), ("isoform_fraction >= 1%", 0.01)]

fig, axes = plt.subplots(1, 2, figsize=(16, 7))
for ax, (label, threshold) in zip(axes, PANELS):
    df = base[base["isoform_fraction"] >= threshold].copy()
    n_tx_per_gene = df.groupby("gene_id")["transcript_id"].transform("count")
    df = df[n_tx_per_gene >= 2].copy()
    n_tx_per_gene = n_tx_per_gene[df.index]
    df["weight"] = 1.0 / n_tx_per_gene

    ax.hist(df["deviation_arctanh"], weights=df["weight"], bins=100, color="steelblue", alpha=0.7)
    ax.axvline(0, color="gray", linestyle="--", linewidth=1)
    ax.set_xlabel("arctanh(tx tau_weighted) - arctanh(gene tau_weighted)")
    ax.set_ylabel("gene-weighted count (tx weighted by 1/n_tx of its gene)")
    ax.set_title(
        f"{label}\nn_tx={len(df)}, n_genes={df['gene_id'].nunique()}, total weight={df['weight'].sum():.0f}"
    )

    print(f"panel={label!r} n_tx={len(df)} n_genes={df['gene_id'].nunique()} "
          f"total_weight={df['weight'].sum():.4f}", file=sys.stderr)

fig.suptitle("gtex_v11_tx_protein_coding: arctanh tau deviation, gene-weighted")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)

print(f"-> {OUT_PNG}", file=sys.stderr)
