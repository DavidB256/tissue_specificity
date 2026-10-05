"""Scatter: one point per multi-transcript gene (n_tx>=2 -- single-tx genes are
trivially weighted_kl=0, uninformative here), x=weighted_kl (TPM-share-weighted
average of KL(tx tissue pmf || gene tissue pmf), bits), y=gene_tau_weighted
(gene's own tau_free_weighted, hierclust0.1/max-agg) -- both from
compute_tx_info_content.py's tx_info_content.csv."""
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IN_CSV = sys.argv[1]
OUT_PNG = sys.argv[2]

df = pd.read_csv(IN_CSV).dropna(subset=["weighted_kl", "gene_tau_weighted"])
df = df[df["n_tx"] >= 2]

fig, ax = plt.subplots(figsize=(8, 7))
ax.scatter(df["weighted_kl"], df["gene_tau_weighted"], s=8, alpha=0.1, color="steelblue")
ax.set_xlabel("weighted KL sum (bits)")
ax.set_ylabel("gene tau_weighted (free_cluster_weights, hierclust0.1, max-agg)")
ax.set_title(f"gtex_v11_tx_protein_coding: weighted KL vs. gene tau_weighted, multi-tx genes only (n={len(df)} genes)")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)

print(f"n={len(df)} -> {OUT_PNG}", file=sys.stderr)
