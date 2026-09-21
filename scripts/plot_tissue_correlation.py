"""Visualize the tissue x tissue correlation matrix (heatmap + dendrogram) using the
precomputed linkage from compute_tissue_correlation.py."""
import sys
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use("Agg")
from plot_utils import greek_label

IN_CORR_CSV = sys.argv[1]
IN_LINKAGE_CSV = sys.argv[2]
OUT_PNG = sys.argv[3]
TITLE_SUFFIX = sys.argv[4] if len(sys.argv) > 4 else "log10 TPM, TPM<1->0"

corr_df = pd.read_csv(IN_CORR_CSV, index_col=0)
Z = pd.read_csv(IN_LINKAGE_CSV).values

g = sns.clustermap(corr_df, row_linkage=Z, col_linkage=Z, cmap="vlag", center=0,
                    figsize=(18, 18), xticklabels=True, yticklabels=True)
g.ax_heatmap.set_xticklabels([greek_label(t.get_text()) for t in g.ax_heatmap.get_xticklabels()], fontsize=5)
g.ax_heatmap.set_yticklabels([greek_label(t.get_text()) for t in g.ax_heatmap.get_yticklabels()], fontsize=5)
g.figure.suptitle(f"tissue correlation ({TITLE_SUFFIX}), n={len(corr_df)} tissues", y=1.02)
g.savefig(OUT_PNG, dpi=300)
