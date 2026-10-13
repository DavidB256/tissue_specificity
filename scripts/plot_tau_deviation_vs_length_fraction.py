"""Test: is (tx_tau_weighted - gene_tau_weighted) driven by tx length
relative to total gene exonic length (length_fraction = tx_length /
gene_exonic_length, union of all the gene's exons)? Signed deviation only.
Panel 1: deviation computed after arctanh-transforming tau (Fisher
z-transform; tau clipped to [0, 1-1e-6] first, arctanh(1) undefined).
Panel 2: deviation on raw tau. All quantified transcripts (tx_total>0) of
multi-tx genes (n_tx>=2), not just "significant" ones."""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

IN_GCT = sys.argv[1]
IN_TX_TAU_CSV = sys.argv[2]
IN_INFO_CONTENT_CSV = sys.argv[3]
IN_LENGTH_FRACTION_CSV = sys.argv[4]
OUT_PNG = sys.argv[5]

EPS = 1e-6


def arctanh_tau(tau):
    return np.arctanh(np.clip(tau, 0.0, 1.0 - EPS))


gct = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(gct.columns[2:])
tx_total = gct[tissues].sum(axis=1)
keep = tx_total > 0

tx = pd.DataFrame({"transcript_id": gct["Name"][keep], "gene_id": gct["Description"][keep]})

tx_tau = pd.read_csv(IN_TX_TAU_CSV).rename(columns={
    "gene": "transcript_id", "tau_free_weighted": "tx_tau_weighted",
})[["transcript_id", "tx_tau_weighted"]]
info = pd.read_csv(IN_INFO_CONTENT_CSV)[["gene_id", "gene_tau_weighted"]]
length_fraction = pd.read_csv(IN_LENGTH_FRACTION_CSV)[["transcript_id", "length_fraction"]]

df = tx.merge(tx_tau, on="transcript_id", how="inner")
df = df.merge(info, on="gene_id", how="inner")
df = df.merge(length_fraction, on="transcript_id", how="inner")
df = df.dropna(subset=["tx_tau_weighted", "gene_tau_weighted", "length_fraction"])
df["deviation"] = df["tx_tau_weighted"] - df["gene_tau_weighted"]
df["deviation_arctanh"] = arctanh_tau(df["tx_tau_weighted"]) - arctanh_tau(df["gene_tau_weighted"])

n_tx_per_gene = df.groupby("gene_id")["transcript_id"].transform("count")
df = df[n_tx_per_gene >= 2]

rho_arctanh, pval_arctanh = spearmanr(df["length_fraction"], df["deviation_arctanh"])
rho_raw, pval_raw = spearmanr(df["length_fraction"], df["deviation"])

fig, axes = plt.subplots(1, 2, figsize=(15, 7))

ax = axes[0]
ax.scatter(df["length_fraction"], df["deviation_arctanh"], s=4, alpha=0.15)
ax.axhline(0, color="gray", linestyle="--", linewidth=1)
ax.set_xlabel("length_fraction (tx_length / gene_exonic_length)")
ax.set_ylabel("arctanh(tx tau_weighted) - arctanh(gene tau_weighted)")
ax.set_title(f"arctanh-transformed deviation\nspearman rho={rho_arctanh:.3f} (p={pval_arctanh:.2g})")

ax = axes[1]
ax.scatter(df["length_fraction"], df["deviation"], s=4, alpha=0.15)
ax.axhline(0, color="gray", linestyle="--", linewidth=1)
ax.set_xlabel("length_fraction (tx_length / gene_exonic_length)")
ax.set_ylabel("tx tau_weighted - gene tau_weighted")
ax.set_title(f"raw deviation\nspearman rho={rho_raw:.3f} (p={pval_raw:.2g})")

fig.suptitle(f"gtex_v11_tx_protein_coding: tau deviation vs. relative tx length (n_tx={len(df)})")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)

print(f"n_tx={len(df)} rho_arctanh={rho_arctanh:.4f} pval_arctanh={pval_arctanh:.4g} "
      f"rho_raw={rho_raw:.4f} pval_raw={pval_raw:.4g} -> {OUT_PNG}", file=sys.stderr)
