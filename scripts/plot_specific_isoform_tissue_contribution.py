"""Bar plot: x = 68 tissues (gtex_v11_tx_protein_coding), y = proportion of
total expression (summed TPM across all transcripts of all genes, that
tissue) contributed by transcripts of "specific_isoform" genes (gene_id in
specific_isoform_genes.csv)."""
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IN_GCT = sys.argv[1]
IN_SPECIFIC_ISOFORM_CSV = sys.argv[2]
OUT_PNG = sys.argv[3]

gct = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(gct.columns[2:])

specific_genes = set(pd.read_csv(IN_SPECIFIC_ISOFORM_CSV)["gene_id"])
is_specific = gct["Description"].isin(specific_genes)

total = gct[tissues].sum(axis=0)
specific_total = gct.loc[is_specific, tissues].sum(axis=0)
proportion = (specific_total / total).reindex(tissues)

fig, ax = plt.subplots(figsize=(16, 7))
ax.bar(tissues, proportion.values, color="steelblue")
ax.set_xticks(range(len(tissues)))
ax.set_xticklabels(tissues, rotation=30, ha="right", fontsize=7)
ax.set_ylabel("proportion of total expression from specific_isoform genes")
ax.set_title(
    f"gtex_v11_tx_protein_coding: tissue expression share of {len(specific_genes)} specific_isoform genes"
)
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)

print(f"n_tissues={len(tissues)} n_specific_genes={len(specific_genes)} "
      f"mean_prop={proportion.mean():.4f} min={proportion.min():.4f} max={proportion.max():.4f} -> {OUT_PNG}",
      file=sys.stderr)
