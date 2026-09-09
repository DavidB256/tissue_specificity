"""Does reducing free_weights (optimize_free_weights.py, one DE param per tissue) down
to free_cluster_weights (optimize_free_cluster_weights.py, one DE param per
hierarchical-clustering cluster at DIST_THRESHOLD=0.1) cost much of the achieved
Spearman rho between weighted_tau(full tissues) and tau(reduced tissues)? For each
scenario, compares the two methods' best-of-N-chains rho (from their chains_summary.csv)
against their optimization dimensionality (n_tissues vs n_clusters).
Args: repeating groups of (label, n_tissues, n_clusters, full_weights_summary_csv,
cluster_weights_summary_csv), then OUT_CSV, OUT_PNG."""
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT_CSV = sys.argv[-2]
OUT_PNG = sys.argv[-1]
GROUPS = sys.argv[1:-2]
assert len(GROUPS) % 5 == 0, "expected repeating groups of 5 args"
SCENARIOS = [GROUPS[i:i + 5] for i in range(0, len(GROUPS), 5)]

rows = []
for label, n_tissues, n_clusters, full_summary_csv, cluster_summary_csv in SCENARIOS:
    n_tissues, n_clusters = int(n_tissues), int(n_clusters)
    full_rho = pd.read_csv(full_summary_csv)["best_rho"].max()
    cluster_rho = pd.read_csv(cluster_summary_csv)["best_rho"].max()
    rows.append({"scenario": label, "method": "free_weights", "n_dims": n_tissues, "best_rho": full_rho})
    rows.append({"scenario": label, "method": "free_cluster_weights", "n_dims": n_clusters, "best_rho": cluster_rho})
    delta = full_rho - cluster_rho
    print(f"{label}: free_weights rho={full_rho:.6f} (d={n_tissues}) vs "
          f"free_cluster_weights rho={cluster_rho:.6f} (d={n_clusters}) -- delta={delta:.6f}", file=sys.stderr)

df = pd.DataFrame(rows)
df.to_csv(OUT_CSV, index=False)

labels = [s[0] for s in SCENARIOS]
full_rhos = df[df["method"] == "free_weights"]["best_rho"].values
cluster_rhos = df[df["method"] == "free_cluster_weights"]["best_rho"].values
n_tissues_list = df[df["method"] == "free_weights"]["n_dims"].values
n_clusters_list = df[df["method"] == "free_cluster_weights"]["n_dims"].values

fig, ax = plt.subplots(figsize=(8, 5))
x = range(len(SCENARIOS))
width = 0.35
ax.bar([i - width / 2 for i in x], full_rhos, width, label="free_weights (per-tissue)", color="tab:blue")
ax.bar([i + width / 2 for i in x], cluster_rhos, width, label="free_cluster_weights (per-cluster)", color="tab:orange")
for i in x:
    ax.text(i - width / 2, full_rhos[i] + 0.002, f"{n_tissues_list[i]}d", ha="center", fontsize=8)
    ax.text(i + width / 2, cluster_rhos[i] + 0.002, f"{n_clusters_list[i]}d", ha="center", fontsize=8)
ax.set_xticks(list(x))
ax.set_xticklabels(labels)
ax.set_ylabel("best Spearman rho (weighted_tau(full) vs tau(reduced))")
ax.set_ylim(min(min(full_rhos), min(cluster_rhos)) - 0.02, 1.0)
ax.legend(loc="lower right", fontsize=8)
ax.set_title("free_weights vs free_cluster_weights: optimization dimensionality vs achieved rho")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)
