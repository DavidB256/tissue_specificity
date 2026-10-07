"""Bar plot, 2 rows: x = 68 tissues (gtex_v11_tx_protein_coding), y = proportion
of total expression (summed TPM across all transcripts of all genes, that
tissue) contributed by transcripts of "specific_isoform" genes (gene_id in
specific_isoform_genes.csv). Bars colored by hierclust0.1 tissue cluster.
Row 1: tissues ordered by cluster, then proportion within cluster.
Row 2: tissues ordered by proportion alone."""
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib as mpl

IN_GCT = sys.argv[1]
IN_SPECIFIC_ISOFORM_CSV = sys.argv[2]
IN_CLUSTERS = sys.argv[3]
OUT_PNG = sys.argv[4]

gct = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(gct.columns[2:])

specific_genes = set(pd.read_csv(IN_SPECIFIC_ISOFORM_CSV)["gene_id"])
is_specific = gct["Description"].isin(specific_genes)

total = gct[tissues].sum(axis=0)
specific_total = gct.loc[is_specific, tissues].sum(axis=0)
proportion = (specific_total / total).reindex(tissues)

clusters = pd.read_csv(IN_CLUSTERS)
tissue_to_cluster = dict(zip(clusters["tissue"], clusters["cluster_id"]))
cluster_label = (clusters.assign(namelen=clusters["tissue"].str.len())
                 .sort_values("namelen")
                 .groupby("cluster_id")["tissue"].first().to_dict())

cluster_size = dict(zip(clusters["cluster_id"], clusters["cluster_size"]))
cluster_ids = sorted(cluster_label, key=lambda c: -cluster_size[c])
palette = mpl.colormaps["tab20"].resampled(len(cluster_ids))
cluster_color = {c: palette(i) for i, c in enumerate(cluster_ids)}

cluster_rank = {c: i for i, c in enumerate(cluster_ids)}
order_by_cluster = sorted(tissues, key=lambda t: (cluster_rank[tissue_to_cluster[t]], proportion[t]))
order_by_proportion = proportion.sort_values(ascending=False).index.tolist()

handles = [plt.Rectangle((0, 0), 1, 1, color=cluster_color[c]) for c in cluster_ids]
legend_labels = [cluster_label[c] for c in cluster_ids]

fig, axes = plt.subplots(2, 1, figsize=(16, 14))
for ax, order, subtitle in [
    (axes[0], order_by_cluster, "ordered by cluster, then proportion"),
    (axes[1], order_by_proportion, "ordered by proportion"),
]:
    bar_colors = [cluster_color[tissue_to_cluster[t]] for t in order]
    ax.bar(order, proportion[order].values, color=bar_colors)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, rotation=30, ha="right", fontsize=7)
    ax.set_ylabel("proportion of total expression from specific_isoform genes")
    ax.set_title(subtitle)
    ax.legend(handles, legend_labels, title="hierclust0.1 cluster", fontsize=6, ncol=2, loc="upper right")

fig.suptitle(f"gtex_v11_tx_protein_coding: tissue expression share of {len(specific_genes)} specific_isoform genes")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)

print(f"n_tissues={len(tissues)} n_specific_genes={len(specific_genes)} "
      f"mean_prop={proportion.mean():.4f} min={proportion.min():.4f} max={proportion.max():.4f} -> {OUT_PNG}",
      file=sys.stderr)
