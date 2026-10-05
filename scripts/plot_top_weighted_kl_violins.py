"""Violin plot: one violin per gene (the top-20 by weighted_kl,
summarize_top_weighted_kl_genes.py's top20_weighted_kl_genes.csv, same
weighted_kl-descending order), values = that gene's per-transcript KL
divergence (KL(tx tissue pmf || gene tissue pmf), bits, from
compute_tx_info_content.py's tx_kl_divergence.pkl). Individual tx KL values
overlaid as jittered points, since some genes have as few as 3 transcripts."""
import sys
import pickle
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IN_CSV = sys.argv[1]
IN_PICKLE = sys.argv[2]
OUT_PNG = sys.argv[3]

top = pd.read_csv(IN_CSV)
with open(IN_PICKLE, "rb") as f:
    kl_dict = pickle.load(f)

labels = []
data = []
for _, row in top.iterrows():
    kls = [v[1] for v in kl_dict[row["gene_id"]].values()]
    data.append(kls)
    labels.append(f"{row['gene_name']}\n(n={row['n_tx']})")

fig, ax = plt.subplots(figsize=(14, 7))
positions = np.arange(1, len(data) + 1)
parts = ax.violinplot(data, positions=positions, showmedians=True, widths=0.8)
for pc in parts["bodies"]:
    pc.set_facecolor("steelblue")
    pc.set_alpha(0.4)

rng = np.random.default_rng(0)
for pos, kls in zip(positions, data):
    jitter = rng.uniform(-0.12, 0.12, size=len(kls))
    ax.scatter(np.full(len(kls), pos) + jitter, kls, s=14, color="steelblue", alpha=0.7, zorder=3)

ax.set_xticks(positions)
ax.set_xticklabels(labels, rotation=30, ha="right")
ax.set_ylabel("per-transcript KL(tx tissue pmf || gene tissue pmf) (bits)")
ax.set_title("gtex_v11_tx_protein_coding: per-tx KL divergence, top 20 genes by weighted_kl")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)

print(f"genes={len(data)} -> {OUT_PNG}", file=sys.stderr)
