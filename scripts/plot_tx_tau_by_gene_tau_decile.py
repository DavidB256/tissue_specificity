"""Violin plots of gene tau_weighted bins ([0.2,0.3), [0.3,0.4), ..., [0.8,0.9),
[0.9,0.95), [0.95,1.0] -- <0.2 dropped, always empty; top decile split at 0.95),
gtex_v11_tx_protein_coding (already gene-level protein_coding filtered in the
base GCT). Each violin: tx tau_weighted for "significant" transcripts -- txs
contributing >=10% of their gene's total expression (isoform_fraction =
tx_total / gene_total TPM, from the base GCT, tx_total>0 only). n annotated
per bin (can be 0)."""
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

SIG_FRACTION = 0.10

gct = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(gct.columns[2:])
tx_total = gct[tissues].sum(axis=1)
keep = tx_total > 0

isoform_fraction_df = pd.DataFrame({
    "transcript_id": gct["Name"][keep],
    "gene_id": gct["Description"][keep],
    "tx_total": tx_total[keep],
})
isoform_fraction_df["gene_total"] = isoform_fraction_df.groupby("gene_id")["tx_total"].transform("sum")
isoform_fraction_df["isoform_fraction"] = isoform_fraction_df["tx_total"] / isoform_fraction_df["gene_total"]

tx_tau = pd.read_csv(IN_TX_TAU_CSV).rename(columns={
    "gene": "transcript_id", "tau_free_weighted": "tx_tau_weighted",
})[["transcript_id", "tx_tau_weighted"]]
info = pd.read_csv(IN_INFO_CONTENT_CSV)[["gene_id", "gene_tau_weighted"]]

df = isoform_fraction_df.merge(tx_tau, on="transcript_id", how="inner")
df = df.merge(info, on="gene_id", how="inner")
df = df.dropna(subset=["isoform_fraction", "tx_tau_weighted", "gene_tau_weighted"])

sig = df[df["isoform_fraction"] >= SIG_FRACTION].copy()

bin_edges = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0 + 1e-9]
n_bins = len(bin_edges) - 1
bin_labels = [f"[{bin_edges[i]:.2f},{bin_edges[i+1]:.2f})" for i in range(n_bins - 1)] + ["[0.95,1.00]"]
sig["bin"] = pd.cut(sig["gene_tau_weighted"], bins=bin_edges, right=False, labels=range(n_bins))

groups = [sig.loc[sig["bin"] == b, "tx_tau_weighted"].dropna().values for b in range(n_bins)]
ns = [len(g) for g in groups]

fig, ax = plt.subplots(figsize=(14, 7))
positions = [i + 1 for i, g in enumerate(groups) if len(g) > 0]
nonempty = [g for g in groups if len(g) > 0]
if nonempty:
    parts = ax.violinplot(nonempty, positions=positions, showmedians=True, widths=0.7)
    for pc in parts["bodies"]:
        pc.set_facecolor("steelblue")
        pc.set_alpha(0.5)

ax.set_xticks(range(1, n_bins + 1))
ax.set_xticklabels([f"{lab}\nn={n}" for lab, n in zip(bin_labels, ns)], fontsize=8)
ax.set_xlabel("gene tau_weighted decile bin")
ax.set_ylabel("tx tau_weighted")
ax.set_title(
    f"gtex_v11_tx_protein_coding: tx tau_weighted by gene tau_weighted decile\n"
    f"(significant tx = isoform_fraction >= {SIG_FRACTION}, n_tx={sum(ns)})"
)
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)

print(f"n_sig_tx={sum(ns)} bin_counts={ns} -> {OUT_PNG}", file=sys.stderr)
