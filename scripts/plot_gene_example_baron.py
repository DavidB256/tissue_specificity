"""4-panel figure using one gene (default INS) as a representative example of
tissue(cell-type)-specificity in the Baron pancreas data, and of why naive
mean-collapse of related cell types (e.g. the 4 islet subtypes -> "oneislets")
understates specificity relative to full-resolution tau, while free_weights
(optimize_free_weights.py, MODE=tissue) largely recovers it.
Panels: (A) per-cell-type TPM at full resolution, (B) tau_unweighted(full) vs
tau_unweighted(oneislets mean-collapse) vs tau_free_weighted(full) for the gene,
(C) genome-wide scatter of tau_unweighted(full) vs tau_unweighted(oneislets),
gene highlighted, (D) per-cell-type free_weights weight, islet subtypes vs rest."""
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

GENE = sys.argv[1] if len(sys.argv) > 1 else "INS"
BARON_DIR = "/ccb/salz4-4/dbass13/expr_data/gene_specificity/baron_protein_coding"
OUT_PNG = sys.argv[2] if len(sys.argv) > 2 else f"{BARON_DIR}/figures/gene_example_{GENE}.png"
ISLET_MEMBERS = ["alpha", "beta", "gamma", "delta"]
COLLAPSED_NAME = "oneislets"

full_gct = pd.read_csv(f"{BARON_DIR}/celltype_tpm_protein_coding.gct.gz", sep="\t", skiprows=2)
full_tau = pd.read_csv(f"{BARON_DIR}/cluster_weights/weighted_tau_simplified_0.1.csv")
collapsed_tau = pd.read_csv(f"{BARON_DIR}/oneislets/weighted_tau.csv")
free_tau = pd.read_csv(f"{BARON_DIR}/free_weights/oneislets/weighted_tau.csv")
free_weights = pd.read_csv(f"{BARON_DIR}/free_weights/oneislets/weights.csv")

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
ax.bar(tpm.index, tpm.values, color=colors_a)
ax.set_yscale("log")
ax.set_ylabel(f"{GENE} TPM (log scale)")
ax.tick_params(axis="x", labelsize=8, rotation=90)
ax.set_title(f"A: {GENE} expression per cell type (full resolution)")
ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color="tab:red", label="islet subtypes"),
                    plt.Rectangle((0, 0), 1, 1, color="tab:gray", label="other")],
          fontsize=8)

# B: tau comparison
ax = axes[0, 1]
labels = [f"tau_unweighted\n(full, n={len(celltypes)})",
          f"tau_unweighted\n({COLLAPSED_NAME}, mean-collapse, n={len(celltypes) - len(ISLET_MEMBERS) + 1})",
          "tau_free_weighted\n(full, weights fit to\noneislets reference)"]
values = [tau_full, tau_collapsed, tau_free]
bar_colors = ["tab:blue", "tab:orange", "tab:green"]
ax.bar(labels, values, color=bar_colors)
ax.set_ylim(0, 1)
ax.set_ylabel("tau")
for i, v in enumerate(values):
    ax.text(i, v + 0.02, f"{v:.3f}", ha="center", fontsize=10)
ax.set_title(f"B: {GENE} specificity: naive collapse loses it, free weights recover it")

# C: genome-wide scatter, tau_full vs tau_collapsed, gene highlighted
ax = axes[1, 0]
merged = full_tau[["gene", "tau_unweighted"]].rename(columns={"tau_unweighted": "tau_full"}).merge(
    collapsed_tau[["gene", "tau_unweighted"]].rename(columns={"tau_unweighted": "tau_collapsed"}), on="gene")
ax.scatter(merged["tau_full"], merged["tau_collapsed"], s=3, alpha=0.15, color="tab:gray", rasterized=True)
ax.plot([0, 1], [0, 1], "--", color="black", linewidth=1)
ax.scatter([tau_full], [tau_collapsed], s=200, color="tab:red", marker="*", zorder=5, label=GENE)
ax.set_xlabel("tau_unweighted (full resolution)")
ax.set_ylabel(f"tau_unweighted ({COLLAPSED_NAME}, mean-collapse)")
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.legend(fontsize=9)
ax.set_title(f"C: genome-wide, {GENE} highlighted (n={len(merged)} genes)")

# D: free_weights per cell type
ax = axes[1, 1]
fw = free_weights.set_index("tissue").loc[celltypes].reset_index()
colors_d = ["tab:red" if ct in ISLET_MEMBERS else "tab:gray" for ct in fw["tissue"]]
ax.bar(fw["tissue"], fw["weight"], color=colors_d)
ax.set_ylim(0, 1)
ax.set_ylabel("free weight")
ax.tick_params(axis="x", labelsize=8, rotation=90)
ax.set_title(f"D: free_weights per cell type (fit to {COLLAPSED_NAME} reference)")

fig.suptitle(f"{GENE} as a representative cell-type-specific gene (Baron pancreas data)", fontsize=13)
fig.tight_layout(rect=(0, 0, 1, 0.97))
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)

print(f"gene={GENE} tau_full={tau_full:.4f} tau_{COLLAPSED_NAME}={tau_collapsed:.4f} "
      f"tau_free={tau_free:.4f} top_celltype={tpm.index[0]} top_tpm={tpm.iloc[0]:.1f}", file=sys.stderr)
print(f"wrote {OUT_PNG}", file=sys.stderr)
