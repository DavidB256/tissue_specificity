"""Shared plotting helpers."""
from scipy.stats import pearsonr, spearmanr

LABEL_DISPLAY_NAMES = {
    "onebrain": "brain tissues combined (56 tissues)",
    "oneislets": "islet cell types combined (α, β, γ, δ, ε)",
}

# Baron pancreas islet subtypes: greek symbol for axis/title brevity. No-op for
# any other tissue/cell-type name (e.g. all GTEx tissues, other Baron types).
CELLTYPE_GREEK = {
    "alpha": "α", "beta": "β", "gamma": "γ", "delta": "δ", "epsilon": "ε",
}


def display_label(label):
    """Human-readable text for a tissue-reduction scenario label (e.g. "onebrain") used
    in figure titles/axis labels."""
    return LABEL_DISPLAY_NAMES.get(label, label)


def greek_label(name):
    """Greek symbol for a Baron islet subtype name, unchanged otherwise."""
    return CELLTYPE_GREEK.get(name, name)


def greek_labels(names):
    return [greek_label(n) for n in names]


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
