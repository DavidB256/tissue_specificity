"""2-panel figure for compare_dendrogram_scenarios.py's output: left, bar chart
of n_clusters (flatclust threshold 0.1) for each of the 8 {pearson,spearman} x
{TPM<1,TPM<10} x {gtex_v11_tx_protein_coding,chess_tx_protein_coding} scenarios;
right, heatmap of the pairwise dendrogram distance matrix (1 - cophenetic
correlation, restricted to common tissues for cross-source pairs)."""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IN_CLUSTER_COUNTS_CSV = sys.argv[1]
IN_DISTANCE_MATRIX_CSV = sys.argv[2]
OUT_PNG = sys.argv[3]

SOURCE_LABEL = {"gtex_v11_tx_protein_coding": "gtex", "chess_tx_protein_coding": "chess"}
SOURCE_COLOR = {"gtex_v11_tx_protein_coding": "#4C72B0", "chess_tx_protein_coding": "#DD8452"}
FLOOR_HATCH = {1.0: "", 10.0: "//"}


def short_name(scenario):
    source, floor, method = scenario.split("__")
    return f"{SOURCE_LABEL[source]}\nf{floor.replace('floor', '')}\n{method}"


counts = pd.read_csv(IN_CLUSTER_COUNTS_CSV)
counts["scenario"] = (counts["source"] + "__floor" + counts["floor_threshold"].astype(int).astype(str)
                       + "__" + counts["corr_method"])
order = ["gtex_v11_tx_protein_coding", "chess_tx_protein_coding"]
counts["source_rank"] = counts["source"].map({s: i for i, s in enumerate(order)})
counts = counts.sort_values(["source_rank", "floor_threshold", "corr_method"])

dist = pd.read_csv(IN_DISTANCE_MATRIX_CSV, index_col=0)
dist = dist.loc[counts["scenario"], counts["scenario"]]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6.5))

x = np.arange(len(counts))
colors = counts["source"].map(SOURCE_COLOR)
hatches = counts["floor_threshold"].map(FLOOR_HATCH)
bars = ax1.bar(x, counts["n_clusters"], color=colors)
for bar, hatch in zip(bars, hatches):
    bar.set_hatch(hatch)
    bar.set_edgecolor("black")
for xi, n in zip(x, counts["n_clusters"]):
    ax1.text(xi, n + 0.5, str(n), ha="center", fontsize=9)
ax1.set_xticks(x)
ax1.set_xticklabels([short_name(s) for s in counts["scenario"]], fontsize=8)
ax1.set_ylabel("n clusters (flatclust threshold 0.1)")
ax1.set_title("tissue cluster count by scenario")
from matplotlib.patches import Patch
legend_elems = [Patch(facecolor=SOURCE_COLOR[s], label=SOURCE_LABEL[s]) for s in order]
legend_elems += [Patch(facecolor="white", edgecolor="black", hatch=h, label=f"TPM<{int(f)}")
                  for f, h in FLOOR_HATCH.items()]
ax1.legend(handles=legend_elems, fontsize=8, loc="upper left")

labels = [short_name(s).replace("\n", " ") for s in counts["scenario"]]
im = ax2.imshow(dist.values, cmap="viridis", vmin=0)
ax2.set_xticks(range(len(labels)))
ax2.set_xticklabels(labels, rotation=90, fontsize=7)
ax2.set_yticks(range(len(labels)))
ax2.set_yticklabels(labels, fontsize=7)
for i in range(len(labels)):
    for j in range(len(labels)):
        v = dist.values[i, j]
        ax2.text(j, i, f"{v:.2f}", ha="center", va="center",
                  color="white" if v < dist.values.max() / 2 else "black", fontsize=6)
ax2.set_title("dendrogram distance (1 - cophenetic r)")
fig.colorbar(im, ax=ax2, fraction=0.046, pad=0.04)

fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)

print(f"-> {OUT_PNG}", file=sys.stderr)
