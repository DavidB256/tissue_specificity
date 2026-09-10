"""Vertical bar chart of free_cluster_weights per-tissue weight (broadcast from the
winning per-cluster DE fit, see optimize_free_cluster_weights.py), sorted descending,
colored by hierarchical-clustering cluster_id (IN_CLUSTERS, DIST_THRESHOLD=0.1) so
tissues sharing one cluster share one color. Generic across datasets (GTEx, Baron).
Weight axis is strictly bounded to [0, 1]."""
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IN_WEIGHTS = sys.argv[1]
IN_CLUSTERS = sys.argv[2]
OUT_PNG = sys.argv[3]

weights = pd.read_csv(IN_WEIGHTS)
clusters = pd.read_csv(IN_CLUSTERS)[["tissue", "cluster_id"]]
df = weights.merge(clusters, on="tissue", how="left").sort_values("weight", ascending=False)

cluster_ids = sorted(df["cluster_id"].unique())
cmap = plt.get_cmap("tab20", max(len(cluster_ids), 1))
color_of_cluster = {c: cmap(i) for i, c in enumerate(cluster_ids)}


def colors_of(d):
    return d["cluster_id"].map(color_of_cluster)


legend_handles = [plt.Rectangle((0, 0), 1, 1, color=color_of_cluster[c], label=f"cluster {c}")
                   for c in cluster_ids]
legend_ncol = max(1, len(cluster_ids) // 10)
title = (f"free_cluster_weights: optimized per-cluster weight, broadcast to tissue "
         f"(n={len(df)} tissues, {len(cluster_ids)} clusters, winner chain)")

fig, ax = plt.subplots(figsize=(max(8, 0.3 * len(df)), 8))
ax.bar(df["tissue"], df["weight"], color=colors_of(df))
ax.set_ylim(0, 1)
ax.set_ylabel("free_cluster weight")
ax.set_xlabel("tissue")
ax.tick_params(axis="x", labelsize=7, rotation=90)
ax.legend(handles=legend_handles, loc="upper right", fontsize=6, ncol=legend_ncol)
ax.set_title(title)
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)

print(f"tissues={len(df)} clusters={len(cluster_ids)}", file=sys.stderr)
