"""tau_unweighted vs tau_clustersize (cluster-size-downweighted, k blocks) scatter with
diagonal stats -- see compute_tau_clustersize.py."""
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plot_utils import annotate_diag_stats

IN_TAU = sys.argv[1]
OUT_PNG = sys.argv[2]
K_LABEL = sys.argv[3]

df = pd.read_csv(IN_TAU)[["tau_unweighted", "tau_clustersize"]].dropna()

fig, ax = plt.subplots(figsize=(6, 6))
ax.scatter(df["tau_unweighted"], df["tau_clustersize"], s=4, alpha=0.15, color="tab:green")
ax.plot([0, 1], [0, 1], "--", color="gray", lw=1)
ax.set_xlabel("tau_unweighted")
ax.set_ylabel(f"tau_clustersize (k={K_LABEL})")
ax.set_title(f"cluster-size downweighting (k={K_LABEL}), n={len(df)} genes")
annotate_diag_stats(ax, df["tau_unweighted"], df["tau_clustersize"])
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)
