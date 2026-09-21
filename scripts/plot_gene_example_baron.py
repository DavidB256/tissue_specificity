"""4-panel figure using one gene (default INS) as a representative example of
tissue(cell-type)-specificity in the Baron pancreas data, and of how well
free_weights/free_cluster_weights (optimize_free_weights.py, MODE=tissue/cluster)
fit to a collapsed "oneislets" reference (mean- or max-collapse of the 5 islet
subtypes) recovers full-resolution tau relative to the naive collapsed tau.
Panels: (A) per-cell-type TPM at full resolution, (B) genome-wide scatter of
tau_unweighted vs tau_free_weighted (full res, weights fit to the oneislets
reference), gene highlighted, (C) genome-wide scatter of tau_unweighted(full)
vs tau_unweighted(oneislets collapse), gene highlighted, (D) per-cell-type free
weight, islet subtypes vs rest."""
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plot_utils import greek_labels

GENE = sys.argv[1] if len(sys.argv) > 1 else "INS"
WEIGHT_MODE = sys.argv[2] if len(sys.argv) > 2 else "free_weights"
assert WEIGHT_MODE in ("free_weights", "free_cluster_weights"), WEIGHT_MODE
AGG = sys.argv[3] if len(sys.argv) > 3 else "mean"
assert AGG in ("mean", "max"), AGG
BARON_DIR = "/ccb/salz4-4/dbass13/expr_data/gene_specificity/baron_protein_coding"
OUT_PNG = sys.argv[4] if len(sys.argv) > 4 else f"{BARON_DIR}/figures/gene_example_{GENE}_{WEIGHT_MODE}_{AGG}.png"
ISLET_MEMBERS = ["alpha", "beta", "gamma", "delta", "epsilon"]
COLLAPSED_NAME = "oneislets"
AGG_SUFFIX = "" if AGG == "mean" else "_max"

full_gct = pd.read_csv(f"{BARON_DIR}/celltype_tpm_protein_coding.gct.gz", sep="\t", skiprows=2)
full_tau = pd.read_csv(f"{BARON_DIR}/cluster_weights/weighted_tau_simplified_0.1.csv")
collapsed_tau = pd.read_csv(f"{BARON_DIR}/{COLLAPSED_NAME}{AGG_SUFFIX}/weighted_tau.csv")
free_tau = pd.read_csv(f"{BARON_DIR}/{WEIGHT_MODE}{AGG_SUFFIX}/{COLLAPSED_NAME}/weighted_tau.csv")
free_weights = pd.read_csv(f"{BARON_DIR}/{WEIGHT_MODE}{AGG_SUFFIX}/{COLLAPSED_NAME}/weights.csv")

celltypes = list(full_gct.columns[2:])
row = full_gct.loc[full_gct["Name"] == GENE].iloc[0]
tpm = row[celltypes].astype(float).sort_values(ascending=False)

tau_full = full_tau.loc[full_tau["gene"] == GENE, "tau_unweighted"].iloc[0]
tau_collapsed = collapsed_tau.loc[collapsed_tau["gene"] == GENE, "tau_unweighted"].iloc[0]
tau_free = free_tau.loc[free_tau["gene"] == GENE, "tau_free_weighted"].iloc[0]

fig, axes = plt.subplots(2, 2, figsize=(13, 10))

# A: per-cell-type TPM, full resolution
ax = axes[0, 0]
colors_a = ["tab:red" if ct in ISLET_MEMBERS else "tab:gray" for ct in tpm.index]
ax.bar(greek_labels(tpm.index), tpm.values, color=colors_a)
ax.set_yscale("log")
ax.set_ylabel(f"{GENE} TPM (log scale)")
ax.tick_params(axis="x", labelsize=8, rotation=90)
ax.set_title(f"A: {GENE} expression per cell type (full resolution)")
ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color="tab:red", label="islet subtypes"),
                    plt.Rectangle((0, 0), 1, 1, color="tab:gray", label="other")],
          fontsize=8)

# B: genome-wide scatter, tau_unweighted vs tau_free_weighted (full res), gene highlighted
ax = axes[0, 1]
ax.scatter(free_tau["tau_unweighted"], free_tau["tau_free_weighted"], s=3, alpha=0.15,
           color="tab:gray", rasterized=True)
ax.plot([0, 1], [0, 1], "--", color="black", linewidth=1)
ax.scatter([tau_full], [tau_free], s=200, color="tab:green", marker="*", zorder=5, label=GENE)
ax.set_xlabel("tau_unweighted (full resolution)")
ax.set_ylabel(f"tau_free_weighted (fit to {COLLAPSED_NAME} {AGG}-collapse)")
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.legend(fontsize=9)
ax.set_title(f"B: {WEIGHT_MODE}, {GENE} highlighted (n={len(free_tau)} genes)")

# C: genome-wide scatter, tau_full vs tau_collapsed, gene highlighted
ax = axes[1, 0]
merged = full_tau[["gene", "tau_unweighted"]].rename(columns={"tau_unweighted": "tau_full"}).merge(
    collapsed_tau[["gene", "tau_unweighted"]].rename(columns={"tau_unweighted": "tau_collapsed"}), on="gene")
ax.scatter(merged["tau_full"], merged["tau_collapsed"], s=3, alpha=0.15, color="tab:gray", rasterized=True)
ax.plot([0, 1], [0, 1], "--", color="black", linewidth=1)
ax.scatter([tau_full], [tau_collapsed], s=200, color="tab:red", marker="*", zorder=5, label=GENE)
ax.set_xlabel("tau_unweighted (full resolution)")
ax.set_ylabel(f"tau_unweighted ({COLLAPSED_NAME}, {AGG}-collapse)")
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.legend(fontsize=9)
ax.set_title(f"C: genome-wide, {GENE} highlighted (n={len(merged)} genes)")

# D: free weight per cell type
ax = axes[1, 1]
fw = free_weights.set_index("tissue").loc[celltypes].reset_index()
colors_d = ["tab:red" if ct in ISLET_MEMBERS else "tab:gray" for ct in fw["tissue"]]
ax.bar(greek_labels(fw["tissue"]), fw["weight"], color=colors_d)
ax.set_ylim(0, 1)
ax.set_ylabel("free weight")
ax.tick_params(axis="x", labelsize=8, rotation=90)
ax.set_title(f"D: {WEIGHT_MODE} per cell type ({COLLAPSED_NAME} {AGG}-collapse ref.)")

fig.suptitle(f"{GENE} as a representative cell-type-specific gene (Baron pancreas data, "
             f"{WEIGHT_MODE}, {AGG}-collapse)", fontsize=13)
fig.tight_layout(rect=(0, 0, 1, 0.97))
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)

print(f"gene={GENE} mode={WEIGHT_MODE} agg={AGG} tau_full={tau_full:.4f} "
      f"tau_{COLLAPSED_NAME}={tau_collapsed:.4f} tau_free={tau_free:.4f} "
      f"top_celltype={tpm.index[0]} top_tpm={tpm.iloc[0]:.1f}", file=sys.stderr)
print(f"wrote {OUT_PNG}", file=sys.stderr)
