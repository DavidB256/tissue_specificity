"""Horizontal bar chart of the free_weights per-tissue weight (winner of N DE chains,
see optimize_free_weights.py / combine_free_weights_chains.py), sorted ascending.
Unlike plot_free_weights.py (GTEx-specific Brain/Kidney/Liver highlighting), tissues
are colored by whether they belong to a multi-member cluster at the clustering
threshold used to build the free_weights reference (cluster_size>1, i.e. the tissue
gets collapsed with at least one other) vs a singleton cluster (kept as its own tissue).
Vertical lines mark the [0.01, 1] optimization bounds (see optimize_free_weights.py)."""
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IN_WEIGHTS = sys.argv[1]
IN_CLUSTERS = sys.argv[2]
OUT_PNG = sys.argv[3]

LO_BOUND, HI_BOUND = 0.01, 1.0

weights = pd.read_csv(IN_WEIGHTS)
clusters = pd.read_csv(IN_CLUSTERS)[["tissue", "cluster_size"]]
df = weights.merge(clusters, on="tissue", how="left").sort_values("weight", ascending=True)

collapsed = df["cluster_size"] > 1
colors = collapsed.map({True: "tab:red", False: "tab:blue"})

fig, ax = plt.subplots(figsize=(8, max(6, 0.3 * len(df))))
ax.barh(df["tissue"], df["weight"], color=colors)
ax.axvline(LO_BOUND, ls=":", color="black", lw=1, label=f"bounds [{LO_BOUND}, {HI_BOUND}]")
ax.axvline(HI_BOUND, ls=":", color="black", lw=1)
ax.set_xlabel("free weight")
ax.set_ylabel("tissue")
ax.tick_params(axis="y", labelsize=7)
ax.legend(handles=[
    plt.Rectangle((0, 0), 1, 1, color="tab:red", label="collapsed (cluster_size>1)"),
    plt.Rectangle((0, 0), 1, 1, color="tab:blue", label="not collapsed (singleton)"),
    ax.lines[0],
], loc="lower right", fontsize=8)
ax.set_title(f"free_weights: optimized per-tissue weight (n={len(df)} tissues, winner chain)")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)
print(f"tissues={len(df)} collapsed={int(collapsed.sum())}", file=sys.stderr)
