"""Vertical bar chart of the diag_reg_inv_corr ("diagonal-regularized inverse
correlation") per-tissue weight (weight_diag_reg_inv_corr, see compute_weighted_tau.py
/ gls_tau.diag_reg_lambda_search), sorted descending. Tissues colored by whether they
belong to a multi-member cluster at the clustering threshold (cluster_size>1) vs a
singleton cluster, same convention as plot_free_weights.py. Weight axis is not fixed
to [0, 1] -- these weights aren't bounded that way."""
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plot_utils import greek_labels

IN_WEIGHTS = sys.argv[1]
IN_CLUSTERS = sys.argv[2]
OUT_PNG = sys.argv[3]

weights = pd.read_csv(IN_WEIGHTS)
clusters = pd.read_csv(IN_CLUSTERS)[["tissue", "cluster_size"]]
df = weights.merge(clusters, on="tissue", how="left").sort_values("weight_diag_reg_inv_corr", ascending=False)


def colors_of(d):
    return (d["cluster_size"] > 1).map({True: "tab:red", False: "tab:blue"})


legend_handles = [
    plt.Rectangle((0, 0), 1, 1, color="tab:red", label="collapsed (cluster_size>1)"),
    plt.Rectangle((0, 0), 1, 1, color="tab:blue", label="not collapsed (singleton)"),
]
title = f"diag_reg_inv_corr: per-tissue weight (n={len(df)} tissues, diag_reg_lambda={df['diag_reg_lambda'].iloc[0]:.4g})"

fig, ax = plt.subplots(figsize=(max(8, 0.3 * len(df)), 6))
ax.bar(greek_labels(df["tissue"]), df["weight_diag_reg_inv_corr"], color=colors_of(df))
ax.set_ylabel("diag_reg_inv_corr weight")
ax.set_xlabel("tissue")
ax.tick_params(axis="x", labelsize=7, rotation=90)
ax.legend(handles=legend_handles, loc="upper right", fontsize=8)
ax.set_title(title)
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)

print(f"tissues={len(df)} collapsed={int((df['cluster_size'] > 1).sum())}", file=sys.stderr)
