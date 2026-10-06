"""Histogram: one point per transcript, value = tx_total / gene_total (share of
the gene's total expression, summed across tissues, that this tx composes).
gtex_v11_tx_protein_coding base GCT (already floored, TPM<1 per tissue -> 0).
Only txs with >=1 TPM in >=1 tissue (tx_total>0) are included."""
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IN_GCT = sys.argv[1]
OUT_PNG = sys.argv[2]

df = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(df.columns[2:])
tx_total = df[tissues].sum(axis=1)

keep = tx_total > 0
tx_total = tx_total[keep]
gene_ids = df["Description"][keep]
gene_total = tx_total.groupby(gene_ids).transform("sum")
share = tx_total / gene_total

fig, ax = plt.subplots(figsize=(8, 6))
ax.hist(share, bins=50, color="steelblue", edgecolor="black")
ax.set_xlabel("tx_total TPM / gene_total TPM (transcript share of gene expression)")
ax.set_ylabel("count (transcripts)")
ax.set_title(f"gtex_v11_tx_protein_coding: transcript share of gene expression (n={len(share)})")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)

print(f"n={len(share)} median={share.median():.4f} mean={share.mean():.4f} -> {OUT_PNG}", file=sys.stderr)
