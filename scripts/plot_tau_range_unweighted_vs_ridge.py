"""Per-gene jackknife tau_range: unweighted vs GLS/ridge-weighted (see
compute_tau_jackknife_pertissue.py / compute_tau_jackknife_pertissue_ridge.py). Points
below y=x mean ridge weighting made tau more stable to a dropped tissue."""
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plot_utils import annotate_diag_stats

IN_UNWEIGHTED = sys.argv[1]
IN_RIDGE = sys.argv[2]
OUT_PNG = sys.argv[3]

u = pd.read_csv(IN_UNWEIGHTED)[["gene", "tau_range"]].rename(columns={"tau_range": "tau_range_unweighted"})
r = pd.read_csv(IN_RIDGE)[["gene", "tau_range"]].rename(columns={"tau_range": "tau_range_ridge"})
df = u.merge(r, on="gene").dropna()

fig, ax = plt.subplots(figsize=(6, 6))
ax.scatter(df["tau_range_unweighted"], df["tau_range_ridge"], s=4, alpha=0.15, color="tab:orange")
hi = max(df["tau_range_unweighted"].max(), df["tau_range_ridge"].max())
ax.plot([0, hi], [0, hi], "--", color="gray", lw=1)
ax.set_xlabel("tau_range (unweighted, jackknife)")
ax.set_ylabel("tau_range (gls_ridge, jackknife)")
ax.set_title(f"jackknife tau_range: unweighted vs ridge, n={len(df)} genes")
annotate_diag_stats(ax, df["tau_range_unweighted"], df["tau_range_ridge"])
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)
