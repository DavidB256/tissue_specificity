"""Poster fig: bar chart of each tissue cluster's weight under
free_cluster_weights_max/hierclust0.1 for gtex_protein_coding (cluster_weights.csv
-- one weight per cluster_id, shared by every gene, fit by the free_cluster DE
optimization). Clusters are hardcoded-named by their member tissues (from
cluster_weights/full/tissue_clusters.csv) for poster brevity -- e.g. cluster 3
(11 non-cerebellar brain regions) -> "Brain", cluster 2 (2 cerebellar regions)
-> "Brain2", cluster 1 (Liver/Hepatocyte/Mixed_Cell) -> "Liver", cluster 14
(Liver_Portal_Tract alone) -> "Liver2". Bars sorted by weight, descending.
dpi=600, poster-sized text."""
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IN_CLUSTER_WEIGHTS_CSV = sys.argv[1]
OUT_PNG = sys.argv[2]

CLUSTER_NAME = {
    1: "Liver", 2: "Brain2", 3: "Brain", 4: "Pituitary", 5: "Spleen",
    6: "Lymphocytes", 7: "Blood", 8: "Pancreas", 9: "Islets", 10: "Heart",
    11: "Muscle", 12: "Stomach", 13: "Intestine", 14: "Liver2", 15: "Skin",
    16: "Kidney", 17: "Mixed", 18: "Glands", 19: "Lung", 20: "Adrenal",
    21: "Fibroblasts", 22: "Testis",
}

AXIS_LABEL_FONTSIZE = 18
TITLE_FONTSIZE = 18
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
ax.set_ylabel("cluster weight (free_cluster_weights_max)", fontsize=AXIS_LABEL_FONTSIZE)
ax.set_ylim(0, 1.08)
ax.set_title("gtex_protein_coding: tissue-cluster weights (hierclust0.1 leg)", fontsize=TITLE_FONTSIZE)

fig.tight_layout()
fig.savefig(OUT_PNG, dpi=600)
plt.close(fig)

print(f"-> {OUT_PNG}", file=sys.stderr)
