"""2-panel scatter comparing tau on the default (full-resolution) tissue set to a
reduced-resolution scenario, both sides computed the same way with no GLS anywhere:
(1) tau_unweighted, default vs reduced -- straight off each resolution's own TPM
values; (2) tau_weighted, default vs reduced -- the free/DE-optimized weighting (see
optimize_free_weights.py for the default side, collapse_free_weighted_tau.py for the
reduced side) applied to each resolution's own TPM values."""
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plot_utils import annotate_diag_stats

IN_DEFAULT_TAU = sys.argv[1]
IN_REDUCED_TAU = sys.argv[2]
OUT_PNG = sys.argv[3]
LABEL = sys.argv[4] if len(sys.argv) > 4 else "reduced"

default = pd.read_csv(IN_DEFAULT_TAU)[["gene", "tau_unweighted", "tau_free_weighted"]]
reduced = pd.read_csv(IN_REDUCED_TAU)[["gene", "tau_unweighted", "tau_free_weighted"]] \
    .rename(columns={"tau_unweighted": "tau_unweighted_reduced", "tau_free_weighted": "tau_free_weighted_reduced"})
df = default.merge(reduced, on="gene").dropna()

PANELS = [
    ("tau_unweighted", "tau_unweighted_reduced", "tau_unweighted"),
    ("tau_free_weighted", "tau_free_weighted_reduced", "tau_weighted (free/DE-optimized)"),
]

fig, axes = plt.subplots(1, 2, figsize=(14, 7))
for ax, (xcol, ycol, title) in zip(axes, PANELS):
    ax.scatter(df[xcol], df[ycol], s=4, alpha=0.2, color="tab:blue")
    ax.plot([0, 1], [0, 1], "--", color="gray", lw=1)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel(f"{xcol} (default)")
    ax.set_ylabel(f"{xcol} ({LABEL})")
    ax.set_title(title, fontsize=11)
    annotate_diag_stats(ax, df[xcol], df[ycol])

fig.suptitle(f"default vs {LABEL}: tau_unweighted and tau_weighted, n={len(df)} genes")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)
print(f"genes={len(df)}", file=sys.stderr)
