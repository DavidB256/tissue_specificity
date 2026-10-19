"""5-panel figure for compare_gene_vs_tx_dendrograms.py's output (pearson
only): top-left, bar chart of n_clusters at flatclust threshold 0.1; top-mid,
bar chart of mean_dist_1_minus_r (mean off-diagonal 1-corr distance -- a
threshold-free, continuous measure of "how similar do tissues look overall",
complementing the cutoff-dependent n_clusters bar); right, n_clusters vs.
flatclust distance threshold for all 9 scenarios (gene-centric source
(GENCODE gene) and the 2 transcript-centric sources (GENCODE tx, CHESS tx) --
all 3 are the same GTEx v11 data, differing only in annotation and gene-vs-tx
aggregation -- each at none/TPM<1/TPM<10 flooring); bottom-left, cophenetic
distance heatmap (height-sensitive); bottom-mid, Robinson-Foulds topology
distance heatmap (height-blind)."""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

IN_CLUSTER_COUNTS_CSV = sys.argv[1]
IN_DISTANCE_MATRIX_CSV = sys.argv[2]
IN_RF_MATRIX_CSV = sys.argv[3]
IN_THRESHOLD_CURVE_CSV = sys.argv[4]
OUT_PNG = sys.argv[5]

SOURCE_ORDER = ["gtex_protein_coding", "gtex_v11_tx_protein_coding", "chess_tx_protein_coding"]
SOURCE_LABEL = {
    "gtex_protein_coding": "GENCODE (gene)",
    "gtex_v11_tx_protein_coding": "GENCODE (tx)",
    "chess_tx_protein_coding": "CHESS (tx)",
}
SOURCE_COLOR = {
    "gtex_protein_coding": "#C44E52",
    "gtex_v11_tx_protein_coding": "#4C72B0",
    "chess_tx_protein_coding": "#55A868",
}
KIND_COLOR = {"gene_centric": "#C44E52", "tx_centric": "#4C72B0"}
FLOOR_LABEL = {0.0: "none", 1.0: "TPM<1", 10.0: "TPM<10"}
FLOOR_HATCH = {0.0: "..", 1.0: "", 10.0: "//"}
FLOOR_LINESTYLE = {0.0: ":", 1.0: "-", 10.0: "--"}
DIST_THRESHOLD = 0.1

counts = pd.read_csv(IN_CLUSTER_COUNTS_CSV)
counts["order_rank"] = counts["source"].map({s: i for i, s in enumerate(SOURCE_ORDER)})
counts = counts.sort_values(["order_rank", "floor_threshold"])


def short_name(scenario):
    source, floor = scenario.split("__floor")
    return f"{SOURCE_LABEL[source]}\n{FLOOR_LABEL[float(floor)]}"


dist = pd.read_csv(IN_DISTANCE_MATRIX_CSV, index_col=0).loc[counts["scenario"], counts["scenario"]]
rf = pd.read_csv(IN_RF_MATRIX_CSV, index_col=0).loc[counts["scenario"], counts["scenario"]]
curves = pd.read_csv(IN_THRESHOLD_CURVE_CSV)

fig, axd = plt.subplot_mosaic(
    [["bar_n", "bar_meanr", "curves"], ["coph", "rf", "curves"]],
    figsize=(21, 11),
)

x = np.arange(len(counts))
colors = counts["kind"].map(KIND_COLOR)
hatches = counts["floor_threshold"].map(FLOOR_HATCH)
legend_elems = [Patch(facecolor=c, edgecolor="black", label=k) for k, c in KIND_COLOR.items()]
legend_elems += [Patch(facecolor="white", edgecolor="black", hatch=h, label=FLOOR_LABEL[f])
                  for f, h in FLOOR_HATCH.items()]


def bar_panel(ax, col, ylabel, title, value_fmt):
    bars = ax.bar(x, counts[col], color=colors, edgecolor="black")
    for bar, hatch in zip(bars, hatches):
        bar.set_hatch(hatch)
    ymax = counts[col].max()
    for xi, v in zip(x, counts[col]):
        ax.text(xi, v + ymax * 0.015, value_fmt(v), ha="center", fontsize=9)
    ax.set_xticks(x)
    ax.set_xticklabels([short_name(s) for s in counts["scenario"]], fontsize=9, rotation=30, ha="right")
    ax.set_ylabel(ylabel)
    ax.set_title(title)


bar_panel(axd["bar_n"], "n_clusters", "n clusters (flatclust threshold 0.1)",
           "tissue cluster count by source x TPM floor\n(pearson only)", lambda v: str(int(v)))
axd["bar_n"].legend(handles=legend_elems, fontsize=8, loc="upper left")

bar_panel(axd["bar_meanr"], "mean_dist_1_minus_r", "mean pairwise distance (1 - pearson r)",
           "mean tissue-pair distance by source x TPM floor\n(pearson only, threshold-free)",
           lambda v: f"{v:.3f}")

ax_curves = axd["curves"]
for name, g in curves.groupby("scenario", sort=False):
    g = g.sort_values("threshold")
    source, floor = name.split("__floor")
    ax_curves.step(g["threshold"], g["n_clusters"], where="post",
                    color=SOURCE_COLOR[source], linestyle=FLOOR_LINESTYLE[float(floor)], linewidth=1.6)
ax_curves.axvline(DIST_THRESHOLD, color="gray", linestyle=(0, (1, 1)), linewidth=1.2)
ax_curves.text(DIST_THRESHOLD, ax_curves.get_ylim()[1] * 0.97, " threshold=0.1", fontsize=8, color="gray")
ax_curves.set_xlabel("flatclust distance threshold")
ax_curves.set_ylabel("n clusters")
ax_curves.set_title("n clusters vs. threshold, all 9 scenarios")
source_legend = [Line2D([0], [0], color=SOURCE_COLOR[s], lw=1.6, label=SOURCE_LABEL[s]) for s in SOURCE_ORDER]
floor_legend = [Line2D([0], [0], color="black", lw=1.6, linestyle=ls, label=FLOOR_LABEL[f])
                for f, ls in FLOOR_LINESTYLE.items()]
leg1 = ax_curves.legend(handles=source_legend, fontsize=8, loc="upper right", title="source")
ax_curves.add_artist(leg1)
ax_curves.legend(handles=floor_legend, fontsize=8, loc="center right", title="floor")

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


plot_heatmap(axd["coph"], dist, "cophenetic distance\n(1 - cophenetic r)")
plot_heatmap(axd["rf"], rf, "Robinson-Foulds distance\n(1 - split jaccard)")

fig.suptitle("gene-centric vs. transcript-centric dendrograms, GENCODE vs CHESS, none/TPM<1/TPM<10 floor (pearson)")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)

print(f"-> {OUT_PNG}", file=sys.stderr)
