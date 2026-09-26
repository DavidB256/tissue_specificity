"""Side-by-side heatmap of a tissue x tissue correlation matrix R and its inverse
R^-1 (np.linalg.inv). Not wired into snakefile -- run manually. Rows/cols reordered
by the precomputed hierarchical clustering leaf order (tissue_correlation_linkage.csv),
no tick labels (too many tissues to read). R uses a truncated "Blues" cmap (matches
plot_tissue_correlation_heatmap_only.py); R^-1 uses a diverging cmap centered at 0,
clipped to the 99th percentile of |off-diagonal entries| since a few tissue pairs with
near-collinear correlation blow up individual inverse entries and would otherwise wash
out the rest of the matrix."""
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

CMAP_LOW_CLIP = 0.15
CMAP_HIGH_CLIP = 0.88
INV_CLIP_PCTL = 99

corr_df = pd.read_csv(IN_CORR_CSV, index_col=0)
Z = pd.read_csv(IN_LINKAGE_CSV)[["cluster1", "cluster2", "distance", "n_obs"]].values
order = dendrogram(Z, no_plot=True, labels=list(corr_df.index))["ivl"]
R = corr_df.loc[order, order].values
Rinv = np.linalg.inv(R)

blues = colormaps["Blues"]
cmap_r = colors.LinearSegmentedColormap.from_list(
    "blues_broad", blues(np.linspace(CMAP_LOW_CLIP, CMAP_HIGH_CLIP, 256)))
off_diag = R[~np.eye(len(R), dtype=bool)]
vmin_r, vmax_r = off_diag.min(), off_diag.max()

off_diag_inv = Rinv[~np.eye(len(Rinv), dtype=bool)]
clip = np.percentile(np.abs(off_diag_inv), INV_CLIP_PCTL)

fig, axes = plt.subplots(1, 2, figsize=(20, 9.5))

im0 = axes[0].imshow(R, cmap=cmap_r, vmin=vmin_r, vmax=vmax_r)
axes[0].set_title("R (tissue correlation)", fontsize=13)
div0 = make_axes_locatable(axes[0])
cax0 = div0.append_axes("right", size="4%", pad=0.08)
fig.colorbar(im0, cax=cax0)

im1 = axes[1].imshow(Rinv, cmap="RdBu_r", vmin=-clip, vmax=clip)
axes[1].set_title(f"R^-1 (clipped to {INV_CLIP_PCTL}th pctl of |off-diag|={clip:.2g})", fontsize=13)
div1 = make_axes_locatable(axes[1])
cax1 = div1.append_axes("right", size="4%", pad=0.08)
fig.colorbar(im1, cax=cax1)

for ax in axes:
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

if TITLE:
    fig.suptitle(TITLE, fontsize=14)
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
plt.close(fig)

eigvals = np.linalg.eigvalsh(R)
print(f"tissues={len(R)} R_offdiag=[{vmin_r:.4f},{vmax_r:.4f}] "
      f"eig_range=[{eigvals.min():.4g},{eigvals.max():.4g}] cond={eigvals.max()/eigvals.min():.4g} "
      f"Rinv_clip={clip:.4g}", file=sys.stderr)
print(f"wrote {OUT_PNG}", file=sys.stderr)
