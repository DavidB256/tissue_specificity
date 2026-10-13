"""Test: is (tx_tau_weighted - gene_tau_weighted) driven by tx length
relative to total gene exonic length (length_fraction = tx_length /
gene_exonic_length, union of all the gene's exons)? Deviation computed after
arctanh-transforming tau (Fisher z-transform; tau scaled by 0.9999 first,
arctanh(1) undefined). Panel 1: significant tx = isoform_fraction >= 10%.
Panel 2: >= 1%. Multi-tx genes only (>=2 tx meeting the threshold)."""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from statsmodels.nonparametric.smoothers_lowess import lowess

LOESS_FRAC = 0.2
SCALE = 0.9999
THRESHOLDS = [0.10, 0.01]

IN_GCT = sys.argv[1]
IN_TX_TAU_CSV = sys.argv[2]
IN_INFO_CONTENT_CSV = sys.argv[3]
IN_LENGTH_FRACTION_CSV = sys.argv[4]
OUT_PNG = sys.argv[5]


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
length_fraction = pd.read_csv(IN_LENGTH_FRACTION_CSV)[["transcript_id", "length_fraction"]]

base = tx.merge(tx_tau, on="transcript_id", how="inner")
base = base.merge(info, on="gene_id", how="inner")
base = base.merge(length_fraction, on="transcript_id", how="inner")
base = base.dropna(subset=["tx_tau_weighted", "gene_tau_weighted", "length_fraction", "isoform_fraction"])
base["deviation_arctanh"] = arctanh_tau(base["tx_tau_weighted"]) - arctanh_tau(base["gene_tau_weighted"])

fig, axes = plt.subplots(1, 2, figsize=(15, 7))
for ax, threshold in zip(axes, THRESHOLDS):
    df = base[base["isoform_fraction"] >= threshold]
    n_tx_per_gene = df.groupby("gene_id")["transcript_id"].transform("count")
    df = df[n_tx_per_gene >= 2]

    rho, pval = spearmanr(df["length_fraction"], df["deviation_arctanh"])
    x_line = np.array([df["length_fraction"].min(), df["length_fraction"].max()])

    ax.scatter(df["length_fraction"], df["deviation_arctanh"], s=4, alpha=0.15)
    ax.axhline(0, color="gray", linestyle="--", linewidth=1)
    slope, intercept = np.polyfit(df["length_fraction"], df["deviation_arctanh"], 1)
    ax.plot(x_line, slope * x_line + intercept, color="red", linewidth=1.5)
    loess_fit = lowess(df["deviation_arctanh"], df["length_fraction"], frac=LOESS_FRAC)
    ax.plot(loess_fit[:, 0], loess_fit[:, 1], color="darkorange", linewidth=1.5)
    ax.set_xlabel("length_fraction (tx_length / gene_exonic_length)")
    ax.set_ylabel("arctanh(tx tau_weighted) - arctanh(gene tau_weighted)")
    ax.set_title(f"significant tx = isoform_fraction >= {threshold:.0%} (n_tx={len(df)})\n"
                 f"spearman rho={rho:.3f} (p={pval:.2g}), OLS slope={slope:.3f}")

    print(f"threshold={threshold} n_tx={len(df)} rho={rho:.4f} pval={pval:.4g} slope={slope:.4f}", file=sys.stderr)

fig.suptitle("gtex_v11_tx_protein_coding: arctanh tau deviation vs. relative tx length")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)

print(f"-> {OUT_PNG}", file=sys.stderr)
