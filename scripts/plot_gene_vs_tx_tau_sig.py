"""Scatter: gene tau_weighted vs. tx tau_weighted, one point per transcript,
for transcripts with isoform_fraction (tx_total / gene_total TPM, tx_total>0
only) > 10%."""
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

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

sig = df[df["isoform_fraction"] > SIG_FRACTION]

rho, pval = spearmanr(sig["gene_tau_weighted"], sig["tx_tau_weighted"])

fig, ax = plt.subplots(figsize=(8, 8))
ax.scatter(sig["gene_tau_weighted"], sig["tx_tau_weighted"], s=4, alpha=0.2)
ax.plot([0, 1], [0, 1], color="gray", linestyle="--", linewidth=1)
ax.set_xlabel("gene tau_weighted")
ax.set_ylabel("tx tau_weighted")
ax.set_title(
    f"gtex_v11_tx_protein_coding: gene vs. tx tau_weighted (isoform_fraction > {SIG_FRACTION})\n"
    f"n_tx={len(sig)}, spearman rho={rho:.3f} (p={pval:.2g})"
)
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)

print(f"n_tx={len(sig)} spearman_rho={rho:.4f} pval={pval:.4g} -> {OUT_PNG}", file=sys.stderr)
