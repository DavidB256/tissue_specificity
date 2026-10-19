"""Hierarchical cluster tree (tissue_correlation_linkage.csv, distance=1-corr, average
linkage) with the flat clusters (cut at DIST_THRESHOLD, see compute_tau_cluster_weighted.py)
highlighted, and (unless SHOW_DIAMETERS=0) each multi-tissue cluster's diameter (max
pairwise distance among its members -- see compute_tau_cluster_weighted.py) drawn as a
horizontal red bar spanning that cluster's leaves. The diameter bar can sit well above
the cluster's own merge height in the tree: average linkage merges on averaged distances,
so a cluster can be admitted by the cut even when some pairs inside it are far more
distant than the cut threshold (chaining) -- this figure makes that visible directly.
Optional trailing args: SHOW_LEGEND (7th, "0" to drop the legend), LW_MULT (8th,
multiplies the tree branch line widths, e.g. 2.0 for a poster render), DPI (9th,
default 300), HEIGHT_MULT (10th, scales figure height, e.g. 0.5 for half as
tall), SHOW_XTICKS (11th, "0" to drop x-axis ticks/leaf labels), SHOW_TITLE
(12th, "0" to drop the suptitle), SHOW_YLABEL (13th, "0" to drop the y-axis
label), LW_ADD (14th, added to the branch line width on top of LW_MULT,
points)."""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram
from plot_utils import display_label, greek_label

IN_LINKAGE_CSV = sys.argv[1]
IN_CLUSTERS_CSV = sys.argv[2]
OUT_PNG = sys.argv[3]
LABEL = sys.argv[4] if len(sys.argv) > 4 else "full"
DIST_THRESHOLD = float(sys.argv[5]) if len(sys.argv) > 5 else 0.1
SHOW_DIAMETERS = sys.argv[6] != "0" if len(sys.argv) > 6 else True
SHOW_LEGEND = sys.argv[7] != "0" if len(sys.argv) > 7 else True
LW_MULT = float(sys.argv[8]) if len(sys.argv) > 8 else 1.0
DPI = int(sys.argv[9]) if len(sys.argv) > 9 else 300
HEIGHT_MULT = float(sys.argv[10]) if len(sys.argv) > 10 else 1.0
SHOW_XTICKS = sys.argv[11] != "0" if len(sys.argv) > 11 else True
SHOW_TITLE = sys.argv[12] != "0" if len(sys.argv) > 12 else True
SHOW_YLABEL = sys.argv[13] != "0" if len(sys.argv) > 13 else True
LW_ADD = float(sys.argv[14]) if len(sys.argv) > 14 else 0.0

clusters = pd.read_csv(IN_CLUSTERS_CSV)
tissues = clusters["tissue"].tolist()
Z = pd.read_csv(IN_LINKAGE_CSV)[["cluster1", "cluster2", "distance", "n_obs"]].values

fig, ax = plt.subplots(figsize=(max(14, 0.22 * len(tissues)), 7 * HEIGHT_MULT))
R = dendrogram(Z, labels=tissues, color_threshold=DIST_THRESHOLD,
               above_threshold_color="lightgray", leaf_rotation=90, leaf_font_size=6, ax=ax)
if LW_MULT != 1.0 or LW_ADD != 0.0:
    for line in ax.get_lines():
        line.set_linewidth(line.get_linewidth() * LW_MULT + LW_ADD)

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

if SHOW_XTICKS:
    ax.set_xticklabels([greek_label(t.get_text()) for t in ax.get_xticklabels()])
else:
    ax.set_xticks([])
ax.axhline(DIST_THRESHOLD, ls="--", color="black", lw=1, label=f"cut = {DIST_THRESHOLD}")
if SHOW_YLABEL:
    ax.set_ylabel("distance (1 - corr)")
if SHOW_LEGEND:
    ax.legend(fontsize=8, loc="upper right")
if SHOW_TITLE:
    fig.suptitle(f"{display_label(LABEL)}: hierarchical tissue clustering, cut={DIST_THRESHOLD}, "
                 f"{clusters['cluster_id'].nunique()} clusters", fontsize=12)
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=DPI)
plt.close(fig)
