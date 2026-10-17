"""3-panel figure for compare_dendrogram_scenarios.py's output: left, bar chart
of n_clusters (flatclust threshold 0.1 -- not a meaningful cut for the UBERON
scenario, see that script's docstring) for each of the 13 scenarios: 12 from
{pearson,spearman} x {none,TPM<1,TPM<10} x {gtex_v11_tx_protein_coding,
chess_tx_protein_coding}, plus 1 from UBERON ontology distance; middle, heatmap
of the pairwise cophenetic distance matrix (height-sensitive); right, heatmap
of the pairwise Robinson-Foulds topology distance matrix (height-blind, pure
tree-shape comparison -- see compare_dendrogram_scenarios.py's docstring for
what distinguishes the two)."""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IN_CLUSTER_COUNTS_CSV = sys.argv[1]
IN_DISTANCE_MATRIX_CSV = sys.argv[2]
IN_RF_MATRIX_CSV = sys.argv[3]
OUT_PNG = sys.argv[4]

SOURCE_LABEL = {"gtex_v11_tx_protein_coding": "gtex", "chess_tx_protein_coding": "chess",
                "gtex_uberon_tree": "uberon"}
SOURCE_COLOR = {"gtex_v11_tx_protein_coding": "#4C72B0", "chess_tx_protein_coding": "#DD8452",
                 "gtex_uberon_tree": "#55A868"}
FLOOR_LABEL = {0.0: "none", 1.0: "TPM<1", 10.0: "TPM<10"}
FLOOR_HATCH = {0.0: "..", 1.0: "", 10.0: "//"}


def short_name(scenario):
    if scenario == "gtex_uberon_tree":
        return "uberon\ntree\n(gtex tissues)"
    source, floor, method = scenario.split("__")
    floor_val = float(floor.replace("floor", ""))
    return f"{SOURCE_LABEL[source]}\n{FLOOR_LABEL[floor_val]}\n{method}"


counts = pd.read_csv(IN_CLUSTER_COUNTS_CSV)
order = ["gtex_v11_tx_protein_coding", "chess_tx_protein_coding", "gtex_uberon_tree"]
counts["source_rank"] = counts["source"].map({s: i for i, s in enumerate(order)})
counts = counts.sort_values(["source_rank", "floor_threshold", "corr_method"])

dist = pd.read_csv(IN_DISTANCE_MATRIX_CSV, index_col=0)
dist = dist.loc[counts["scenario"], counts["scenario"]]
rf = pd.read_csv(IN_RF_MATRIX_CSV, index_col=0)
rf = rf.loc[counts["scenario"], counts["scenario"]]

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(22, 6.5))

x = np.arange(len(counts))
colors = counts["source"].map(SOURCE_COLOR)
hatches = counts["floor_threshold"].map(FLOOR_HATCH).fillna("xx")
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
legend_elems += [Patch(facecolor="white", edgecolor="black", hatch=h, label=FLOOR_LABEL[f])
                  for f, h in FLOOR_HATCH.items()]
ax1.legend(handles=legend_elems, fontsize=8, loc="upper left")

labels = [short_name(s).replace("\n", " ") for s in counts["scenario"]]


def plot_heatmap(ax, mat, title):
    im = ax.imshow(mat.values, cmap="viridis", vmin=0)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=90, fontsize=7)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=7)
    vmax = np.nanmax(mat.values)
    for i in range(len(labels)):
        for j in range(len(labels)):
            v = mat.values[i, j]
            text = "nan" if np.isnan(v) else f"{v:.2f}"
            color = "black" if np.isnan(v) or v < vmax / 2 else "white"
            ax.text(j, i, text, ha="center", va="center", color=color, fontsize=6)
    ax.set_title(title)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)


plot_heatmap(ax2, dist, "cophenetic distance (1 - cophenetic r)\nheight-sensitive")
plot_heatmap(ax3, rf, "Robinson-Foulds distance (1 - split jaccard)\nheight-blind, topology only")

fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)

print(f"-> {OUT_PNG}", file=sys.stderr)
