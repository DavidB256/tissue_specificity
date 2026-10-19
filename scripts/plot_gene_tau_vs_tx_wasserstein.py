"""Scatter: gene tau_weighted (free_cluster_weights, hierclust0.1 leg, max-agg
-- from tx_info_content.csv's gene_tau_weighted, the tx-resolution weight
vector applied to the gene's summed-isoform tissue profile) vs mean
Wasserstein-1 distance between a gene's expressed transcripts and the gene
itself (from compute_gene_tx_wasserstein.py). Single-isoform genes
(n_tx==1, mean distance trivially 0) plotted separately from multi-isoform
genes to show where isoform-level tissue-distribution diversity actually
occurs."""
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr

IN_TX_INFO_CONTENT_CSV = sys.argv[1]
IN_WASSERSTEIN_CSV = sys.argv[2]
OUT_PNG = sys.argv[3]

info = pd.read_csv(IN_TX_INFO_CONTENT_CSV)[["gene_id", "gene_tau_weighted"]]
wd = pd.read_csv(IN_WASSERSTEIN_CSV)
df = info.merge(wd, on="gene_id").dropna()

single = df[df["n_tx"] == 1]
multi = df[df["n_tx"] > 1]

fig, ax = plt.subplots(figsize=(8, 7))
ax.scatter(single["gene_tau_weighted"], single["mean_wasserstein_tx_to_gene"],
           s=4, alpha=0.15, color="tab:blue", label=f"single-isoform (n={len(single)})")
ax.scatter(multi["gene_tau_weighted"], multi["mean_wasserstein_tx_to_gene"],
           s=4, alpha=0.25, color="tab:red", label=f"multi-isoform (n={len(multi)})")
ax.set_xlabel("gene tau_weighted")
ax.set_ylabel("mean Wasserstein-1 dist. (tx PMF vs gene PMF)")
ax.set_title("gtex_v11_tx_protein_coding: gene tau_weighted vs mean tx-to-gene distance")
ax.legend(loc="upper left", fontsize=8, markerscale=3)

r, _ = pearsonr(df["gene_tau_weighted"], df["mean_wasserstein_tx_to_gene"])
rho, _ = spearmanr(df["gene_tau_weighted"], df["mean_wasserstein_tx_to_gene"])
r_multi, _ = pearsonr(multi["gene_tau_weighted"], multi["mean_wasserstein_tx_to_gene"])
ax.text(0.97, 0.03,
        f"all genes: pearson r={r:.3f}, spearman ρ={rho:.3f}\nmulti-isoform only: pearson r={r_multi:.3f}",
        transform=ax.transAxes, va="bottom", ha="right", fontsize=8,
        bbox=dict(facecolor="white", alpha=0.7, edgecolor="none"))

fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)

print(f"-> {OUT_PNG}", file=sys.stderr)
