"""3-panel figure for compare_gene_vs_tx_dendrograms.py's output (pearson
only, flat-cut 0.1): left, bar chart of n_clusters for the gene-centric
source (gtex_protein_coding) and the 2 transcript-centric sources
(gtex_v11_tx_protein_coding, chess_tx_protein_coding), each at both TPM<1 and
TPM<10 flooring; middle, cophenetic distance heatmap (height-sensitive);
right, Robinson-Foulds topology distance heatmap (height-blind)."""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

IN_CLUSTER_COUNTS_CSV = sys.argv[1]
IN_DISTANCE_MATRIX_CSV = sys.argv[2]
IN_RF_MATRIX_CSV = sys.argv[3]
OUT_PNG = sys.argv[4]

SOURCE_LABEL = {
    "gtex_protein_coding": "gtex (gene)",
    "gtex_v11_tx_protein_coding": "gtex (tx)",
    "chess_tx_protein_coding": "chess (tx)",
}
KIND_COLOR = {"gene_centric": "#C44E52", "tx_centric": "#4C72B0"}
FLOOR_HATCH = {1.0: "", 10.0: "//"}

counts = pd.read_csv(IN_CLUSTER_COUNTS_CSV)
order = ["gtex_protein_coding", "gtex_v11_tx_protein_coding", "chess_tx_protein_coding"]
counts["order_rank"] = counts["source"].map({s: i for i, s in enumerate(order)})
counts = counts.sort_values(["order_rank", "floor_threshold"])


def short_name(scenario):
    source, floor = scenario.split("__floor")
    return f"{SOURCE_LABEL[source]}\nTPM<{floor}"


dist = pd.read_csv(IN_DISTANCE_MATRIX_CSV, index_col=0).loc[counts["scenario"], counts["scenario"]]
rf = pd.read_csv(IN_RF_MATRIX_CSV, index_col=0).loc[counts["scenario"], counts["scenario"]]

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6.5))

x = np.arange(len(counts))
colors = counts["kind"].map(KIND_COLOR)
hatches = counts["floor_threshold"].map(FLOOR_HATCH)
bars = ax1.bar(x, counts["n_clusters"], color=colors, edgecolor="black")
for bar, hatch in zip(bars, hatches):
    bar.set_hatch(hatch)
for xi, n in zip(x, counts["n_clusters"]):
    ax1.text(xi, n + 0.3, str(n), ha="center", fontsize=10)
ax1.set_xticks(x)
ax1.set_xticklabels([short_name(s) for s in counts["scenario"]], fontsize=9)
ax1.set_ylabel("n clusters (flatclust threshold 0.1)")
ax1.set_title("tissue cluster count by source x TPM floor\n(pearson only)")
legend_elems = [Patch(facecolor=c, edgecolor="black", label=k) for k, c in KIND_COLOR.items()]
legend_elems += [Patch(facecolor="white", edgecolor="black", hatch=h, label=f"TPM<{int(f)}")
                  for f, h in FLOOR_HATCH.items()]
ax1.legend(handles=legend_elems, fontsize=9, loc="upper left")

labels = [short_name(s).replace("\n", " ") for s in counts["scenario"]]


def plot_heatmap(ax, mat, title):
    im = ax.imshow(mat.values, cmap="viridis", vmin=0)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=8)
    finite = mat.values[np.isfinite(mat.values)]
    vmax = finite.max() if len(finite) else 1.0
    for i in range(len(labels)):
        for j in range(len(labels)):
            v = mat.values[i, j]
            text = "n/a" if np.isnan(v) else f"{v:.2f}"
            color = "black" if np.isnan(v) or v < vmax / 2 else "white"
            ax.text(j, i, text, ha="center", va="center", color=color, fontsize=8)
    ax.set_title(title)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)


plot_heatmap(ax2, dist, "cophenetic distance\n(1 - cophenetic r)")
plot_heatmap(ax3, rf, "Robinson-Foulds distance\n(1 - split jaccard)")

fig.suptitle("gene-centric vs. transcript-centric dendrograms, TPM<1 vs TPM<10 floor (pearson)")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)

print(f"-> {OUT_PNG}", file=sys.stderr)
