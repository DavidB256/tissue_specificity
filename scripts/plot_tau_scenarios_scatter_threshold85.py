"""Single-threshold (0.85) variant of plot_tau_scenarios_scatter_thresholds.py's
2-panel scatter, with the same resizing/font/key adjustments made to poster_figs/
tau_free_cluster_vs_full_onebrain.png (see
plot_poster_tau_free_cluster_vs_full_onebrain.py): wider figure (19x9.5in,
"plenty of inches"), panel titles "tau unweighted"/"tau weighted" at title
fontsize 15, axis labels at fontsize 15, pearson r dropped entirely from the
stats key (spearman rho/below line/above line only), the stats key and the
category legend at a matched font size (13), and no overall suptitle.
Categorizes each gene by whether tau crosses 0.85 going from
tau_unweighted(full) to tau_free_weighted(full) -- 3 mutually exclusive
categories (above->below, below->above, none) instead of the original's 5
(2 thresholds: 0.5 and 0.8). dpi set via optional 5th arg (default 300);
optional 6th arg FONT_BUMP (default 0) adds extra pt to every text element
(key/axis-label/title/tick-label) -- used to size the poster_figs copy up
further (+4pt) without changing the default gene_specificity output; optional
7th arg SHOW_TITLE ("0"/"1", default "1") -- pass "0" to drop the per-panel
titles entirely (for a poster fig with no title at all)."""
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from plot_utils import display_label

IN_FREE_TAU = sys.argv[1]
IN_REF_TAU = sys.argv[2]
OUT_PNG = sys.argv[3]
LABEL = sys.argv[4] if len(sys.argv) > 4 else "reference"
DPI = int(sys.argv[5]) if len(sys.argv) > 5 else 300
FONT_BUMP = int(sys.argv[6]) if len(sys.argv) > 6 else 0
SHOW_TITLE = (sys.argv[7] if len(sys.argv) > 7 else "1") != "0"

THRESHOLD = 0.85
KEY_FONTSIZE = 13 + FONT_BUMP
AXIS_LABEL_FONTSIZE = 15 + FONT_BUMP
TITLE_FONTSIZE = 15 + FONT_BUMP
TICK_LABEL_FONTSIZE = 10 + FONT_BUMP

free = pd.read_csv(IN_FREE_TAU)[["gene", "tau_unweighted", "tau_free_weighted"]]
ref = pd.read_csv(IN_REF_TAU)[["gene", "tau_unweighted"]].rename(columns={"tau_unweighted": "tau_ref"})
df = free.merge(ref, on="gene").dropna()

u, w = df["tau_unweighted"], df["tau_free_weighted"]


def categorize(u, w):
    if u > THRESHOLD and w < THRESHOLD:
        return "above 0.85 to below"
    if u < THRESHOLD and w > THRESHOLD:
        return "below 0.85 to above"
    return "no threshold crossed"


df["category"] = [categorize(ui, wi) for ui, wi in zip(u, w)]

CATEGORY_STYLE = [
    ("above 0.85 to below", "tab:red"),
    ("below 0.85 to above", "tab:purple"),
    ("no threshold crossed", "tab:gray"),
]

PANELS = [
    ("tau_unweighted", "tau unweighted", "tau (full)"),
    ("tau_free_weighted", "tau weighted", "tau_weighted (full, free/DE-optimized)"),
]


def annotate_stats(ax, x, y):
    rho, _ = spearmanr(x, y)
    frac_below = (x > y).mean()
    frac_above = (x < y).mean()
    lines = [f"spearman ρ={rho:.3f}", f"below line={frac_below:.1%}", f"above line={frac_above:.1%}"]
    ax.text(0.03, 0.97, "\n".join(lines), transform=ax.transAxes, va="top", ha="left",
            fontsize=KEY_FONTSIZE, bbox=dict(facecolor="white", alpha=0.7, edgecolor="none"))


fig, axes = plt.subplots(1, 2, figsize=(19, 9.5))
for ax, (xcol, title, xlabel) in zip(axes, PANELS):
    for cat, color in CATEGORY_STYLE:
        sub = df[df["category"] == cat]
        alpha = 0.15 if cat == "no threshold crossed" else 0.6
        size = 6 if cat == "no threshold crossed" else 15  # 1.5x of original (4, 10)
        zorder = 2 if cat == "no threshold crossed" else 3
        ax.scatter(sub[xcol], sub["tau_ref"], s=size, alpha=alpha, color=color,
                   label=f"{cat} (n={len(sub)})", zorder=zorder, rasterized=True)
    ax.plot([0, 1], [0, 1], "--", color="gray", lw=1, zorder=1)
    ax.axvline(THRESHOLD, color="lightgray", lw=0.8, zorder=0)
    ax.axhline(THRESHOLD, color="lightgray", lw=0.8, zorder=0)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel(xlabel, fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel(f"tau_unweighted ({display_label(LABEL)})", fontsize=AXIS_LABEL_FONTSIZE)
    if SHOW_TITLE:
        ax.set_title(title, fontsize=TITLE_FONTSIZE)
    ax.tick_params(labelsize=TICK_LABEL_FONTSIZE)
    annotate_stats(ax, df[xcol], df["tau_ref"])

axes[1].legend(fontsize=KEY_FONTSIZE, loc="lower right")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=DPI)
plt.close(fig)

print(f"genes={len(df)}", file=sys.stderr)
for cat, _ in CATEGORY_STYLE:
    print(f"  {cat}: n={(df['category'] == cat).sum()}", file=sys.stderr)
print(f"wrote {OUT_PNG}", file=sys.stderr)
