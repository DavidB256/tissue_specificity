"""Vertical bar chart of the free_weights per-tissue weight (winner of N DE chains,
see optimize_free_weights.py / combine_free_weights_chains.py), sorted descending.
Tissues colored by whether they belong to a multi-member cluster at the clustering
threshold used to build the free_weights reference (cluster_size>1, i.e. the tissue
gets collapsed with at least one other) vs a singleton cluster (kept as its own
tissue). Weight axis is strictly bounded to [0, 1]."""
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IN_WEIGHTS = sys.argv[1]
IN_CLUSTERS = sys.argv[2]
OUT_PNG = sys.argv[3]

weights = pd.read_csv(IN_WEIGHTS)
clusters = pd.read_csv(IN_CLUSTERS)[["tissue", "cluster_size"]]
df = weights.merge(clusters, on="tissue", how="left").sort_values("weight", ascending=False)


def colors_of(d):
    return (d["cluster_size"] > 1).map({True: "tab:red", False: "tab:blue"})


legend_handles = [
    plt.Rectangle((0, 0), 1, 1, color="tab:red", label="collapsed (cluster_size>1)"),
    plt.Rectangle((0, 0), 1, 1, color="tab:blue", label="not collapsed (singleton)"),
]
title = f"free_weights: optimized per-tissue weight (n={len(df)} tissues, winner chain)"

fig, ax = plt.subplots(figsize=(max(8, 0.3 * len(df)), 6))
ax.bar(df["tissue"], df["weight"], color=colors_of(df))
ax.set_ylim(0, 1)
ax.set_ylabel("free weight")
ax.set_xlabel("tissue")
ax.tick_params(axis="x", labelsize=7, rotation=90)
ax.legend(handles=legend_handles, loc="upper right", fontsize=8)
ax.set_title(title)
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)

print(f"tissues={len(df)} collapsed={int((df['cluster_size'] > 1).sum())}", file=sys.stderr)
