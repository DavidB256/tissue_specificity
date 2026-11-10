"""3-panel scatter (2 on top, 1 centered below) comparing a weighted-tau fit (default: free_weights, see
optimize_free_weights.py; also reused for diag_reg_inv_corr, see
compute_weighted_tau.py) to
its reduced-resolution reference, both against tau_unweighted computed directly from
the reference's own (collapsed) TPM values -- no GLS anywhere in the reference: (1)
tau in the two scenarios -- tau_unweighted on the full (full-resolution) tissue set vs
tau_unweighted on the reference; (2) tau_weighted in the two scenarios -- the weighted
tau column (full) vs the same reference tau_unweighted. tau_weighted is only ever
computed at the full resolution; the reduced side has no separate "weighted" tau of
its own. The bottom panel is tau_unweighted (full) vs tau_weighted (full): how far the
fit moves tau away from the unweighted values."""
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plot_utils import annotate_diag_stats, display_label

IN_FREE_TAU = sys.argv[1]
IN_REF_TAU = sys.argv[2]
OUT_PNG = sys.argv[3]
LABEL = sys.argv[4] if len(sys.argv) > 4 else "reference"
WEIGHTED_COL = sys.argv[5] if len(sys.argv) > 5 else "tau_free_weighted"
METHOD_LABEL = sys.argv[6] if len(sys.argv) > 6 else "free/DE-optimized"

free = pd.read_csv(IN_FREE_TAU)[["gene", "tau_unweighted", WEIGHTED_COL]]
ref = pd.read_csv(IN_REF_TAU)[["gene", "tau_unweighted"]].rename(columns={"tau_unweighted": "tau_ref"})
df = free.merge(ref, on="gene").dropna()

PANELS = [
    ("tau_unweighted", "tau", "tau (full)"),
    (WEIGHTED_COL, "tau_weighted", f"tau_weighted (full, {METHOD_LABEL})"),
]

fig = plt.figure(figsize=(14, 14))
gs = fig.add_gridspec(2, 4)
axes = [fig.add_subplot(gs[0, 0:2]), fig.add_subplot(gs[0, 2:4])]
ax_bottom = fig.add_subplot(gs[1, 1:3])
for ax, (xcol, title, xlabel) in zip(axes, PANELS):
    ax.scatter(df[xcol], df["tau_ref"], s=4, alpha=0.2, color="tab:blue")
    ax.plot([0, 1], [0, 1], "--", color="gray", lw=1)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(f"tau_unweighted ({display_label(LABEL)})")
    ax.set_title(title, fontsize=11)
    annotate_diag_stats(ax, df[xcol], df["tau_ref"])

ax_bottom.scatter(df["tau_unweighted"], df[WEIGHTED_COL], s=4, alpha=0.2, color="tab:blue")
ax_bottom.plot([0, 1], [0, 1], "--", color="gray", lw=1)
ax_bottom.set_xlim(0, 1)
ax_bottom.set_ylim(0, 1)
ax_bottom.set_xlabel("tau_unweighted (full)")
ax_bottom.set_ylabel(f"tau_weighted (full, {METHOD_LABEL})")
ax_bottom.set_title("tau_unweighted vs tau_weighted (both full)", fontsize=11)
annotate_diag_stats(ax_bottom, df["tau_unweighted"], df[WEIGHTED_COL])

fig.suptitle(f"full vs {display_label(LABEL)}: tau and tau_weighted, n={len(df)} genes")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)
print(f"genes={len(df)}", file=sys.stderr)
