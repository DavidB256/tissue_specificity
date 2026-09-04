"""2-panel scatter comparing a free_weights fit (see optimize_free_weights.py) to its
reduced-resolution reference: (1) tau in the two scenarios -- tau_unweighted on the
default (full-resolution) tissue set vs tau_unweighted on the collapsed reference set;
(2) tau_weighted in the two scenarios -- tau_free_weighted (default, DE-optimized
against the reference) vs tau_weighted_gls_ridge (the reference's own GLS-ridge-weighted
tau, computed across its collapsed tissues)."""
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plot_utils import annotate_diag_stats

IN_FREE_TAU = sys.argv[1]
IN_REF_TAU = sys.argv[2]
OUT_PNG = sys.argv[3]
LABEL = sys.argv[4] if len(sys.argv) > 4 else "reference"

free = pd.read_csv(IN_FREE_TAU)[["gene", "tau_unweighted", "tau_free_weighted"]]
ref = pd.read_csv(IN_REF_TAU)[["gene", "tau_unweighted", "tau_weighted_gls_ridge"]] \
    .rename(columns={"tau_unweighted": "tau_ref", "tau_weighted_gls_ridge": "tau_weighted_ref"})
df = free.merge(ref, on="gene").dropna()

PANELS = [
    ("tau_unweighted", "tau_ref", "tau", "tau (default)", f"tau ({LABEL})"),
    ("tau_free_weighted", "tau_weighted_ref", "tau_weighted",
     "tau_weighted (default, free/DE-optimized)", f"tau_weighted ({LABEL}, GLS-ridge)"),
]

fig, axes = plt.subplots(1, 2, figsize=(14, 7))
for ax, (xcol, ycol, title, xlabel, ylabel) in zip(axes, PANELS):
    ax.scatter(df[xcol], df[ycol], s=4, alpha=0.2, color="tab:blue")
    ax.plot([0, 1], [0, 1], "--", color="gray", lw=1)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=11)
    annotate_diag_stats(ax, df[xcol], df[ycol])

fig.suptitle(f"default vs {LABEL}: tau and tau_weighted, n={len(df)} genes")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)
print(f"genes={len(df)}", file=sys.stderr)
