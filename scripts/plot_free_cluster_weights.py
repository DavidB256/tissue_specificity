"""Bar chart of free_cluster_weights per-tissue weight (broadcast from the winning
per-cluster DE fit, see optimize_free_cluster_weights.py), colored by
hierarchical-clustering cluster_id (IN_CLUSTERS, DIST_THRESHOLD=0.1) so tissues
sharing one cluster share one color. Generic across datasets (GTEx, Baron). Weight
axis is strictly bounded to [0, 1]. Writes both a horizontal (tissues on y-axis,
sorted ascending) and a vertical (tissues on x-axis, sorted descending) version."""
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IN_WEIGHTS = sys.argv[1]
IN_CLUSTERS = sys.argv[2]
OUT_PNG_H = sys.argv[3]
OUT_PNG_V = sys.argv[4]

weights = pd.read_csv(IN_WEIGHTS)
clusters = pd.read_csv(IN_CLUSTERS)[["tissue", "cluster_id"]]
df = weights.merge(clusters, on="tissue", how="left")

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

df_h = df.sort_values("weight", ascending=True)
fig, ax = plt.subplots(figsize=(8, max(6, 0.3 * len(df))))
ax.barh(df_h["tissue"], df_h["weight"], color=colors_of(df_h))
ax.set_xlim(0, 1)
ax.set_xlabel("free_cluster weight")
ax.set_ylabel("tissue")
ax.tick_params(axis="y", labelsize=7)
ax.legend(handles=legend_handles, loc="lower right", fontsize=6, ncol=legend_ncol)
ax.set_title(title)
fig.tight_layout()
fig.savefig(OUT_PNG_H, dpi=300)
plt.close(fig)

df_v = df.sort_values("weight", ascending=False)
fig, ax = plt.subplots(figsize=(max(8, 0.3 * len(df)), 8))
ax.bar(df_v["tissue"], df_v["weight"], color=colors_of(df_v))
ax.set_ylim(0, 1)
ax.set_ylabel("free_cluster weight")
ax.set_xlabel("tissue")
ax.tick_params(axis="x", labelsize=7, rotation=90)
ax.legend(handles=legend_handles, loc="upper right", fontsize=6, ncol=legend_ncol)
ax.set_title(title)
fig.tight_layout()
fig.savefig(OUT_PNG_V, dpi=300)
plt.close(fig)

print(f"tissues={len(df)} clusters={len(cluster_ids)}", file=sys.stderr)
