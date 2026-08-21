"""Shared plotting helpers."""
import numpy as np
from scipy.stats import pearsonr, spearmanr


def annotate_diag_stats(ax, x, y, loc=(0.03, 0.97), fontsize=9):
    """Annotate a y=x-diagonal comparison scatterplot (same-scale x/y, e.g. two tau
    variants) with pearson r, spearman rho, and the fraction of points below the line
    (x > y) and above the line (x < y). Points exactly on the line (x == y) are excluded
    from both fractions."""
    r, _ = pearsonr(x, y)
    rho, _ = spearmanr(x, y)
    frac_below = (x > y).mean()
    frac_above = (x < y).mean()
    ax.text(loc[0], loc[1],
            f"pearson r={r:.3f}\nspearman ρ={rho:.3f}\nbelow line={frac_below:.1%}\nabove line={frac_above:.1%}",
            transform=ax.transAxes, va="top", ha="left", fontsize=fontsize,
            bbox=dict(facecolor="white", alpha=0.7, edgecolor="none"))
    return r, rho, frac_below, frac_above


def bootstrap_scatter_panel(fig, ax, x, y, c, title, xlabel, ylabel,
                             floor=1e-5, clabel="tau_unweighted (mean across resamples)"):
    """Log-log diagonal-comparison scatter (left-censored at `floor`), colored by `c`,
    with y=x reference line and annotate_diag_stats. Used for both the bootstrap
    variance and coefficient-of-variation figures."""
    x = np.maximum(x, floor)
    y = np.maximum(y, floor)
    sc = ax.scatter(x, y, c=c, cmap="viridis", s=4, alpha=0.1)
    fig.colorbar(sc, ax=ax, label=clabel)
    lo, hi = min(x.min(), y.min()), max(x.max(), y.max())
    ax.plot([lo, hi], [lo, hi], "--", color="gray", lw=1)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=10)
    return annotate_diag_stats(ax, x, y)
