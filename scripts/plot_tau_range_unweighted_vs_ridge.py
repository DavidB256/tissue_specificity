"""Per-gene jackknife tau_range: unweighted vs GLS/ridge-weighted (see
compute_tau_jackknife_pertissue.py / compute_tau_jackknife_pertissue_ridge.py). Points
below y=x mean ridge weighting made tau more stable to a dropped tissue. Colored by
all-tissue (canonical, no tissue dropped) tau_unweighted from weighted_tau.csv."""
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plot_utils import annotate_diag_stats

IN_UNWEIGHTED = sys.argv[1]
IN_RIDGE = sys.argv[2]
IN_TAU = sys.argv[3]
OUT_PNG = sys.argv[4]
MIN_RANGE = float(sys.argv[5]) if len(sys.argv) > 5 else None

u = pd.read_csv(IN_UNWEIGHTED)[["gene", "tau_range"]].rename(columns={"tau_range": "tau_range_unweighted"})
r = pd.read_csv(IN_RIDGE)[["gene", "tau_range"]].rename(columns={"tau_range": "tau_range_ridge"})
tau = pd.read_csv(IN_TAU)[["gene", "tau_unweighted"]]
df = u.merge(r, on="gene").merge(tau, on="gene").dropna()
if MIN_RANGE is not None:
    df = df[(df["tau_range_unweighted"] >= MIN_RANGE) & (df["tau_range_ridge"] >= MIN_RANGE)]

fig, ax = plt.subplots(figsize=(6, 6))
sc = ax.scatter(df["tau_range_unweighted"], df["tau_range_ridge"], c=df["tau_unweighted"],
                 cmap="viridis", s=4, alpha=0.1)
fig.colorbar(sc, ax=ax, label="tau_unweighted (all-tissue)")
hi = max(df["tau_range_unweighted"].max(), df["tau_range_ridge"].max())
ax.plot([0, hi], [0, hi], "--", color="gray", lw=1)
ax.set_xlabel("tau_range (unweighted, jackknife)")
ax.set_ylabel("tau_range (gls_ridge, jackknife)")
title_suffix = f", tau_range>={MIN_RANGE} both" if MIN_RANGE is not None else ""
ax.set_title(f"jackknife tau_range: unweighted vs ridge, n={len(df)} genes{title_suffix}")
annotate_diag_stats(ax, df["tau_range_unweighted"], df["tau_range_ridge"])
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)
