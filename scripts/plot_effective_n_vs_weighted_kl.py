"""Scatter: one point per gene, x=effective_n (2^H, "effective number of
transcripts in play"), y=weighted_kl (TPM-share-weighted average of
KL(tx tissue pmf || gene tissue pmf), bits) -- both from
compute_tx_info_content.py's tx_info_content.csv."""
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IN_CSV = sys.argv[1]
OUT_PNG = sys.argv[2]

df = pd.read_csv(IN_CSV)

fig, ax = plt.subplots(figsize=(8, 7))
ax.scatter(df["effective_n"], df["weighted_kl"], s=8, alpha=0.2, color="steelblue")
ax.set_xlabel(r"$n^{tx}_{eff} = 2^H$ (effective number of transcripts)")
ax.set_ylabel("weighted KL sum (bits)")
ax.set_title(f"gtex_v11_tx_protein_coding: effective n transcripts vs. weighted KL (n={len(df)} genes)")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)

print(f"n={len(df)} -> {OUT_PNG}", file=sys.stderr)
