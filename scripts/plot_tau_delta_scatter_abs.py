"""Absolute-value version of plot_tau_delta_scatter.py: scatter |delta_tau| vs
|delta_tau_w| for each tissue-reduction scenario (onebrain, simplified_0.1). Points
below the y=x line have |delta_tau_w| < |delta_tau| -- free-weighting shrank the gene's
tau shift from collapsing tissue resolution, regardless of the shift's sign."""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plot_utils import annotate_diag_stats

IN_ONEBRAIN_FREE_TAU = sys.argv[1]
IN_ONEBRAIN_REF_TAU = sys.argv[2]
IN_SIMPLIFIED_FREE_TAU = sys.argv[3]
IN_SIMPLIFIED_REF_TAU = sys.argv[4]
OUT_PNG = sys.argv[5]

SCENARIOS = [
    ("onebrain", IN_ONEBRAIN_FREE_TAU, IN_ONEBRAIN_REF_TAU),
    ("simplified_0.1", IN_SIMPLIFIED_FREE_TAU, IN_SIMPLIFIED_REF_TAU),
]


def load_abs_deltas(free_tau_csv, ref_tau_csv):
    d = pd.read_csv(free_tau_csv)[["gene", "tau_unweighted", "tau_free_weighted"]]
    r = pd.read_csv(ref_tau_csv)[["gene", "tau_unweighted"]].rename(columns={"tau_unweighted": "tau_ref"})
    df = d.merge(r, on="gene").dropna()
    delta_tau = (df["tau_unweighted"] - df["tau_ref"]).abs()
    delta_tau_w = (df["tau_free_weighted"] - df["tau_ref"]).abs()
    return delta_tau, delta_tau_w


fig, axes = plt.subplots(1, len(SCENARIOS), figsize=(7 * len(SCENARIOS), 7))
for ax, (label, free_tau_csv, ref_tau_csv) in zip(axes, SCENARIOS):
    delta_tau, delta_tau_w = load_abs_deltas(free_tau_csv, ref_tau_csv)

    lim = np.percentile(np.concatenate([delta_tau, delta_tau_w]), 99)
    ax.scatter(delta_tau, delta_tau_w, s=4, alpha=0.1, color="tab:blue")
    ax.plot([0, lim], [0, lim], "--", color="gray", lw=1)
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_xlabel(r"$|\Delta\tau|$")
    ax.set_ylabel(r"$|\Delta\tau_w|$")
    ax.set_title(f"default vs {label}, n={len(delta_tau)} genes", fontsize=10)
    annotate_diag_stats(ax, delta_tau, delta_tau_w)

fig.suptitle(r"$|\Delta\tau|$ vs $|\Delta\tau_w|$: magnitude of tau change from collapsing tissue resolution, unweighted vs free-weighted")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)
