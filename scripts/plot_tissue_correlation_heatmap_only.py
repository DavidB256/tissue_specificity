"""Bare heatmap variant of plot_tissue_correlation.py: rows/cols reordered by the
precomputed hierarchical clustering leaf order (tissue_correlation_linkage.csv) but
with no dendrogram tree drawn, no tick labels. Colormap is a truncated "Blues" (skips
the near-white bottom slice of the colormap so low-but-still-positive correlations
don't wash out toward white), with vmin/vmax set to the actual off-diagonal data range
to maximize contrast between cells. Colorbar is a separate axis, exactly as tall as the
heatmap."""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colormaps, colors
from scipy.cluster.hierarchy import dendrogram
from mpl_toolkits.axes_grid1 import make_axes_locatable

IN_CORR_CSV = sys.argv[1]
IN_LINKAGE_CSV = sys.argv[2]
OUT_PNG = sys.argv[3]
TITLE = sys.argv[4] if len(sys.argv) > 4 else None

CMAP_LOW_CLIP = 0.15  # skip the near-white bottom slice of "Blues"
CMAP_HIGH_CLIP = 0.88  # skip the near-black top slice of "Blues"

corr_df = pd.read_csv(IN_CORR_CSV, index_col=0)
Z = pd.read_csv(IN_LINKAGE_CSV)[["cluster1", "cluster2", "distance", "n_obs"]].values
order = dendrogram(Z, no_plot=True, labels=list(corr_df.index))["ivl"]
mat = corr_df.loc[order, order].values

off_diag = mat[~np.eye(len(mat), dtype=bool)]
vmin, vmax = off_diag.min(), off_diag.max()

blues = colormaps["Blues"]
cmap = colors.LinearSegmentedColormap.from_list(
    "blues_broad", blues(np.linspace(CMAP_LOW_CLIP, CMAP_HIGH_CLIP, 256)))

fig, ax = plt.subplots(figsize=(12, 12))
im = ax.imshow(mat, cmap=cmap, vmin=vmin, vmax=vmax)
ax.set_xticks([])
ax.set_yticks([])
for spine in ax.spines.values():
    spine.set_visible(False)
if TITLE:
    ax.set_title(TITLE, fontsize=14)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="4%", pad=0.08)
fig.colorbar(im, cax=cax)

fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"tissues={len(mat)} vmin={vmin:.4f} vmax={vmax:.4f}", file=sys.stderr)
