"""Poster fig: bar chart of each tissue cluster's weight under
free_cluster_weights_max/hierclust0.1 for gtex_protein_coding_noRedundancy (cluster_weights.csv
-- one weight per cluster_id, shared by every gene, fit by the free_cluster DE
optimization). Clusters are hardcoded-named by their member tissues (from
cluster_weights/tissue_clusters_hierclust0.1.csv) for poster brevity -- e.g. cluster 5
(11 non-cerebellar brain regions) -> "Brain2", cluster 4 (2 cerebellar regions)
-> "Brain" (names swapped so Brain2 sits right of Brain). Bars sorted by weight, descending.
dpi=600, poster-sized text."""
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IN_CLUSTER_WEIGHTS_CSV = sys.argv[1]
OUT_PNG = sys.argv[2]

CLUSTER_NAME = {
    1: "Spleen", 2: "Blood", 3: "Lymphocytes", 4: "Brain", 5: "Brain2",
    6: "Pituitary", 7: "Heart", 8: "Muscle", 9: "Skin", 10: "Kidney",
    11: "Intestine", 12: "Glands", 13: "Mixed", 14: "Lung", 15: "Adrenal",
    16: "Fibroblasts", 17: "Pancreas", 18: "Liver", 19: "Testis",
}

AXIS_LABEL_FONTSIZE = 18
TICK_LABEL_FONTSIZE = 16
BAR_LABEL_FONTSIZE = 14

df = pd.read_csv(IN_CLUSTER_WEIGHTS_CSV)
df["cluster_name"] = df["cluster_id"].map(CLUSTER_NAME)
df = df.sort_values("weight", ascending=False).reset_index(drop=True)

fig, ax = plt.subplots(figsize=(16, 8))
x = range(len(df))
bars = ax.bar(x, df["weight"], color="tab:blue", edgecolor="black")
for xi, v in zip(x, df["weight"]):
    ax.text(xi, v + 0.015, f"{v:.2f}", ha="center", fontsize=BAR_LABEL_FONTSIZE)
ax.set_xticks(x)
ax.set_xticklabels(df["cluster_name"], fontsize=TICK_LABEL_FONTSIZE, rotation=30, ha="right",
                    rotation_mode="anchor")
ax.tick_params(axis="y", labelsize=TICK_LABEL_FONTSIZE)
ax.set_xlabel("Tissue cluster", fontsize=AXIS_LABEL_FONTSIZE)
ax.set_ylabel("weight of tissues in cluster", fontsize=AXIS_LABEL_FONTSIZE)
ax.set_ylim(0, 1)

fig.tight_layout()
fig.savefig(OUT_PNG, dpi=600)
plt.close(fig)

print(f"-> {OUT_PNG}", file=sys.stderr)
