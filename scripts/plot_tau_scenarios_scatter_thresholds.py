"""Variant of plot_tau_scenarios_scatter.py that additionally color-labels each gene by
whether tau crosses the 0.5 or 0.8 threshold when going from tau_unweighted(full) to
tau_free_weighted(full) (i.e. comparing the two panels' x-values to each other, not to
the reference/y-axis). A gene crossing both thresholds in the same direction is counted
under the more extreme (0.8) category only, so the 5 categories are mutually exclusive.
Adds thin light-gray reference lines at 0.5 and 0.8 on both axes of both panels."""
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

free = pd.read_csv(IN_FREE_TAU)[["gene", "tau_unweighted", "tau_free_weighted"]]
ref = pd.read_csv(IN_REF_TAU)[["gene", "tau_unweighted"]].rename(columns={"tau_unweighted": "tau_ref"})
df = free.merge(ref, on="gene").dropna()

u, w = df["tau_unweighted"], df["tau_free_weighted"]


def categorize(u, w):
    if u > 0.8 and w < 0.8:
        return "above 0.8 to below"
    if u < 0.8 and w > 0.8:
        return "below 0.8 to above"
    if u > 0.5 and w < 0.5:
        return "above 0.5 to below"
    if u < 0.5 and w > 0.5:
        return "below 0.5 to above"
    return "no threshold crossed"


df["category"] = [categorize(ui, wi) for ui, wi in zip(u, w)]

CATEGORY_STYLE = [
    ("above 0.5 to below", "tab:orange"),
    ("below 0.5 to above", "tab:blue"),
    ("above 0.8 to below", "tab:red"),
    ("below 0.8 to above", "tab:purple"),
    ("no threshold crossed", "tab:gray"),
]

PANELS = [
    ("tau_unweighted", "tau", "tau (full)"),
    ("tau_free_weighted", "tau_weighted", "tau_weighted (full, free/DE-optimized)"),
]

fig, axes = plt.subplots(1, 2, figsize=(15, 7.5))
for ax, (xcol, title, xlabel) in zip(axes, PANELS):
    for cat, color in CATEGORY_STYLE:
        sub = df[df["category"] == cat]
        alpha = 0.15 if cat == "no threshold crossed" else 0.6
        size = 4 if cat == "no threshold crossed" else 10
        zorder = 2 if cat == "no threshold crossed" else 3
        ax.scatter(sub[xcol], sub["tau_ref"], s=size, alpha=alpha, color=color,
                   label=f"{cat} (n={len(sub)})", zorder=zorder, rasterized=True)
    ax.plot([0, 1], [0, 1], "--", color="gray", lw=1, zorder=1)
    for thresh in (0.5, 0.8):
        ax.axvline(thresh, color="lightgray", lw=0.8, zorder=0)
        ax.axhline(thresh, color="lightgray", lw=0.8, zorder=0)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(f"tau_unweighted ({display_label(LABEL)})")
    ax.set_title(title, fontsize=11)
    annotate_diag_stats(ax, df[xcol], df["tau_ref"])

axes[1].legend(fontsize=8, loc="lower right")
fig.suptitle(f"full vs {display_label(LABEL)}: tau and tau_weighted, colored by "
             f"tau_unweighted(full)->tau_weighted(full) threshold crossing, n={len(df)} genes")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)

print(f"genes={len(df)}", file=sys.stderr)
for cat, _ in CATEGORY_STYLE:
    print(f"  {cat}: n={(df['category'] == cat).sum()}", file=sys.stderr)
print(f"wrote {OUT_PNG}", file=sys.stderr)
