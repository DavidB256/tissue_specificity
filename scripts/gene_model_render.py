"""Simple IGV-style gene model renderer: gene span + every transcript's
exon/CDS structure, no track stacking logic or sequence/coverage display --
just coordinates and shapes. Used by plot_extreme_tau_deviation_screenshots.py."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

HIGHLIGHT_COLOR = "crimson"
OTHER_COLOR = "steelblue"
EXON_HEIGHT = 0.3
CDS_HEIGHT = 0.6
OUTLINE_COLOR = "black"
OUTLINE_WIDTH = 1.0
VLINE_COLOR = "grey"


def _axis_scale(span_bp):
    if span_bp >= 1_000_000:
        return 1_000_000.0, "Mb"
    return 1_000.0, "kb"


def render_gene_model(gene, transcripts, order, highlight_tx_ids, vlines, out_png, subtitle=None):
    """gene: dict with gene_id, gene_name, chrom, start, end, strand.
    transcripts: dict transcript_id -> {"strand": str, "exons": [(s,e),...],
    "cds": [(s,e),...]}. order: transcript_ids top-to-bottom (caller sorts,
    e.g. by abundance). highlight_tx_ids: set of transcript_ids drawn in a
    different color. vlines: list of (x_bp, label) grey vertical lines (e.g.
    TSS/5'SS of the highlighted tx). out_png: output path (saved at
    dpi=300)."""
    n_tx = len(order)
    scale, unit = _axis_scale(gene["end"] - gene["start"] + 1)

    fig_height = max(2.5, 0.6 * n_tx + 1.5)
    fig, ax = plt.subplots(figsize=(12, fig_height))

    gene_y = n_tx + 0.8
    ax.plot([gene["start"] / scale, gene["end"] / scale], [gene_y, gene_y],
            color="black", linewidth=2, zorder=1)
    n_arrows = max(1, int((gene["end"] - gene["start"]) / scale * 4))
    for j in range(1, n_arrows + 1):
        x = gene["start"] + (gene["end"] - gene["start"]) * j / (n_arrows + 1)
        marker = ">" if gene["strand"] == "+" else "<"
        ax.plot(x / scale, gene_y, marker=marker, color="black", markersize=5, zorder=2)
    ax.text((gene["start"] + gene["end"]) / 2 / scale, gene_y + 0.3,
            f"{gene['gene_name']} ({gene['gene_id']}) {gene['chrom']}:{gene['start']}-{gene['end']} "
            f"({gene['strand']})",
            ha="center", va="bottom", fontsize=9, fontweight="bold")

    for i, tx_id in enumerate(order):
        y = n_tx - i
        tx = transcripts[tx_id]
        color = HIGHLIGHT_COLOR if tx_id in highlight_tx_ids else OTHER_COLOR
        exon_start = min(s for s, _ in tx["exons"])
        exon_end = max(e for _, e in tx["exons"])

        ax.plot([exon_start / scale, exon_end / scale], [y, y], color=color, linewidth=1, zorder=1)

        for s, e in tx["exons"]:
            ax.add_patch(Rectangle((s / scale, y - EXON_HEIGHT / 2), (e - s) / scale, EXON_HEIGHT,
                                    facecolor=color, edgecolor=OUTLINE_COLOR, linewidth=OUTLINE_WIDTH,
                                    alpha=0.5, zorder=2))
        for s, e in tx.get("cds", []):
            ax.add_patch(Rectangle((s / scale, y - CDS_HEIGHT / 2), (e - s) / scale, CDS_HEIGHT,
                                    facecolor=color, edgecolor=OUTLINE_COLOR, linewidth=OUTLINE_WIDTH,
                                    alpha=1.0, zorder=3))

        ax.text(gene["start"] / scale - (gene["end"] - gene["start"]) / scale * 0.01, y, tx_id,
                ha="right", va="center", fontsize=7,
                fontweight="bold" if tx_id in highlight_tx_ids else "normal")

    for x_bp, label in vlines:
        ax.axvline(x_bp / scale, color=VLINE_COLOR, linestyle="--", linewidth=1, alpha=0.8, zorder=4)
        ax.text(x_bp / scale, gene_y + 0.1, label, ha="center", va="bottom", fontsize=6,
                color=VLINE_COLOR, rotation=90)

    ax.set_xlim(gene["start"] / scale - (gene["end"] - gene["start"]) / scale * 0.25,
                gene["end"] / scale + (gene["end"] - gene["start"]) / scale * 0.02)
    ax.set_ylim(0, gene_y + 1)
    ax.set_yticks([])
    ax.set_xlabel(f"{gene['chrom']} position ({unit})")
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)

    if subtitle:
        ax.set_title(subtitle, fontsize=9)

    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    plt.close(fig)
