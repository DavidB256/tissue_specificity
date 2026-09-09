"""Hierarchical cluster tree (tissue_correlation_linkage.csv, distance=1-corr, average
linkage) with the flat clusters (cut at DIST_THRESHOLD, see compute_tau_cluster_weighted.py)
highlighted, and (unless SHOW_DIAMETERS=0) each multi-tissue cluster's diameter (max
pairwise distance among its members -- see compute_tau_cluster_weighted.py) drawn as a
horizontal red bar spanning that cluster's leaves. The diameter bar can sit well above
the cluster's own merge height in the tree: average linkage merges on averaged distances,
so a cluster can be admitted by the cut even when some pairs inside it are far more
distant than the cut threshold (chaining) -- this figure makes that visible directly."""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram
from plot_utils import display_label

IN_LINKAGE_CSV = sys.argv[1]
IN_CLUSTERS_CSV = sys.argv[2]
OUT_PNG = sys.argv[3]
LABEL = sys.argv[4] if len(sys.argv) > 4 else "full"
DIST_THRESHOLD = float(sys.argv[5]) if len(sys.argv) > 5 else 0.1
SHOW_DIAMETERS = sys.argv[6] != "0" if len(sys.argv) > 6 else True

clusters = pd.read_csv(IN_CLUSTERS_CSV)
tissues = clusters["tissue"].tolist()
Z = pd.read_csv(IN_LINKAGE_CSV)[["cluster1", "cluster2", "distance", "n_obs"]].values

fig, ax = plt.subplots(figsize=(max(14, 0.22 * len(tissues)), 7))
R = dendrogram(Z, labels=tissues, color_threshold=DIST_THRESHOLD,
               above_threshold_color="lightgray", leaf_rotation=90, leaf_font_size=6, ax=ax)

leaf_x = {t: 5 + 10 * i for i, t in enumerate(R["ivl"])}

if SHOW_DIAMETERS:
    for cid, grp in clusters.groupby("cluster_id"):
        if len(grp) < 2:
            continue
        xs = [leaf_x[t] for t in grp["tissue"]]
        diam = grp["diameter"].iloc[0]
        ax.hlines(diam, min(xs), max(xs), color="red", lw=1.5, zorder=5)
        ax.text((min(xs) + max(xs)) / 2, diam, f"{diam:.3f}", fontsize=6, ha="center",
                va="bottom", color="red", zorder=5)
    ax.plot([], [], color="red", lw=1.5, label="cluster diameter (max pairwise distance)")

ax.axhline(DIST_THRESHOLD, ls="--", color="black", lw=1, label=f"cut = {DIST_THRESHOLD}")
ax.set_ylabel("distance (1 - corr)")
ax.legend(fontsize=8, loc="upper right")
fig.suptitle(f"{display_label(LABEL)}: hierarchical tissue clustering, cut={DIST_THRESHOLD}, "
             f"{clusters['cluster_id'].nunique()} clusters", fontsize=12)
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)
