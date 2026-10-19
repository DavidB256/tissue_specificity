"""Standalone 2-panel fig (just the two bar charts from plot_gene_vs_tx_
dendrograms.py), wider for x-label room, "none" floor excluded (TPM<1 filt. /
TPM<10 filt. only). Left: n_clusters at flatclust threshold 0.1 (read from
compare_gene_vs_tx_dendrograms.py's cluster_counts CSV). Right: mean pairwise
tissue distance (1-pearson r), with 1 SE error bars from a nonparametric
bootstrap -- 200 reps, each resampling tissues with replacement and
recomputing the mean off-diagonal distance of the resulting tissue x tissue
submatrix (repeated-tissue pairs contribute distance 0, as is standard for
this kind of unit-resampling bootstrap); SE = std of the 200 resampled means."""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import MultipleLocator
from gls_tau import tissue_corr_expressed_only, floor_tpm

IN_GTEX_PC_GCT = sys.argv[1]
IN_GTEX_TX_PC_GCT = sys.argv[2]
IN_CHESS_TX_PC_GCT = sys.argv[3]
IN_CLUSTER_COUNTS_CSV = sys.argv[4]
OUT_PNG = sys.argv[5]

N_BOOT = 200
RNG = np.random.default_rng(0)
FLOORS = [1.0, 10.0]
METHOD = "pearson"
SOURCE_ORDER = ["gtex_protein_coding", "gtex_v11_tx_protein_coding", "chess_tx_protein_coding"]
SOURCE_LABEL = {
    "gtex_protein_coding": "GENCODE (gene)",
    "gtex_v11_tx_protein_coding": "GENCODE (tx)",
    "chess_tx_protein_coding": "CHESS (tx)",
}
SOURCE_PATH = {
    "gtex_protein_coding": IN_GTEX_PC_GCT,
    "gtex_v11_tx_protein_coding": IN_GTEX_TX_PC_GCT,
    "chess_tx_protein_coding": IN_CHESS_TX_PC_GCT,
}
KIND = {"gtex_protein_coding": "gene_centric", "gtex_v11_tx_protein_coding": "tx_centric",
        "chess_tx_protein_coding": "tx_centric"}
KIND_LABEL = {"gene_centric": "Gene-level", "tx_centric": "Tx-level"}
KIND_COLOR = {"gene_centric": "#C44E52", "tx_centric": "#4C72B0"}
FLOOR_LABEL = {1.0: "TPM<1 filt.", 10.0: "TPM<10 filt."}
FLOOR_HATCH = {1.0: "", 10.0: "//"}


def short_name(source, floor):
    return f"{SOURCE_LABEL[source]}\n{FLOOR_LABEL[floor]}"


def bootstrap_mean_dist(dist, n_boot=N_BOOT, rng=RNG):
    n = dist.shape[0]
    iu = np.triu_indices(n, k=1)
    boots = np.empty(n_boot)
    for b in range(n_boot):
        idx = rng.integers(0, n, size=n)
        boots[b] = dist[np.ix_(idx, idx)][iu].mean()
    return boots.std(ddof=1)


counts = pd.read_csv(IN_CLUSTER_COUNTS_CSV)
counts = counts[counts["floor_threshold"].isin(FLOORS)].copy()
counts["order_rank"] = counts["source"].map({s: i for i, s in enumerate(SOURCE_ORDER)})
counts = counts.sort_values(["order_rank", "floor_threshold"]).reset_index(drop=True)

mean_dist_rows = []
for source in SOURCE_ORDER:
    df = pd.read_csv(SOURCE_PATH[source], sep="\t", skiprows=2)
    tissues = list(df.columns[2:])
    raw_mat = df[tissues].to_numpy(dtype=float)
    for floor in FLOORS:
        mat = floor_tpm(raw_mat, floor=floor)
        R = tissue_corr_expressed_only(mat, method=METHOD)
        dist = 1 - R
        np.fill_diagonal(dist, 0.0)
        dist = (dist + dist.T) / 2
        mean_dist = dist[np.triu_indices_from(dist, k=1)].mean()
        se = bootstrap_mean_dist(dist)
        mean_dist_rows.append({"source": source, "floor_threshold": floor,
                                "mean_dist": mean_dist, "boot_se": se})
        print(f"{source} floor={floor}: mean_dist={mean_dist:.4f} boot_se={se:.4f}", file=sys.stderr)

mean_dist_df = pd.DataFrame(mean_dist_rows)
mean_dist_df["order_rank"] = mean_dist_df["source"].map({s: i for i, s in enumerate(SOURCE_ORDER)})
mean_dist_df = mean_dist_df.sort_values(["order_rank", "floor_threshold"]).reset_index(drop=True)

labels = [short_name(s, f) for s, f in zip(counts["source"], counts["floor_threshold"])]
x = np.arange(len(counts))
colors = counts["source"].map(KIND).map(KIND_COLOR)
hatches = counts["floor_threshold"].map(FLOOR_HATCH)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6.5))

bars1 = ax1.bar(x, counts["n_clusters"], color=colors, edgecolor="black")
for bar, hatch in zip(bars1, hatches):
    bar.set_hatch(hatch)
for xi, n in zip(x, counts["n_clusters"]):
    ax1.text(xi, n + 0.3, str(int(n)), ha="center", fontsize=10)
ax1.set_xticks(x)
ax1.set_xticklabels(labels, fontsize=10, rotation=30, ha="right")
ax1.set_ylabel("num. of clusters (fcluster threshold 0.1)")
legend_elems = [Patch(facecolor=c, edgecolor="black", label=KIND_LABEL[k]) for k, c in KIND_COLOR.items()]
legend_elems += [Patch(facecolor="white", edgecolor="black", hatch=h, label=FLOOR_LABEL[f])
                  for f, h in FLOOR_HATCH.items()]
ax1.legend(handles=legend_elems, fontsize=9, loc="upper left")

bars2 = ax2.bar(x, mean_dist_df["mean_dist"], yerr=mean_dist_df["boot_se"], capsize=4,
                 color=colors, edgecolor="black", error_kw={"linewidth": 1.2})
for bar, hatch in zip(bars2, hatches):
    bar.set_hatch(hatch)
ax2.set_xticks(x)
ax2.set_xticklabels(labels, fontsize=10, rotation=30, ha="right")
ax2.set_ylabel("mean pairwise distance (1 - r)")
ax2.yaxis.set_major_locator(MultipleLocator(0.1))

fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)

print(f"-> {OUT_PNG}", file=sys.stderr)
