"""Poster fig, 2 panels vs. TPM threshold T (x: 1..30; retained = median TPM >= T in at
least one tissue): left, proportion of all protein-coding genes / protein-coding
transcripts in GENCODE v47 that are retained; right, mean pairwise tissue distance
(1 - r) with 95% bootstrap band (precomputed by plot_mean_tissue_dist_vs_tpm_floor.py),
gene vs. transcript. Both with 95% bootstrap bands (resampling gene clusters -- the
gene itself for genes, all of a gene's transcripts for transcripts -- with replacement
from the full GENCODE v47 protein-coding universe, N_BOOT=200). Colors: matplotlib defaults (gene blue, transcript orange). dpi=600.
Usage: ... IN_GENE_GCT IN_TX_GCT IN_MEAN_DIST_CSV IN_GENCODE_GTF OUT_PNG GENE_LABEL TX_LABEL"""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IN_GENE_GCT, IN_TX_GCT, IN_DIST_CSV, IN_GTF, OUT_PNG, GENE_SRC, TX_SRC = sys.argv[1:8]

AXIS_LABEL_FONTSIZE = 19
TICK_LABEL_FONTSIZE = 14
LEGEND_FONTSIZE = 17
COLORS = {"gene": "C0", "transcript": "C1"}
THRESHOLDS = np.arange(1, 31)
N_BOOT = 200
LINE_WIDTH = 3.0   # 2x matplotlib default (1.5)
MARKER_SIZE = 12   # 2x matplotlib default (6)
rng = np.random.default_rng(0)

n_gene_total = n_tx_total = 0
with open(IN_GTF) as f:
    for line in f:
        if line[0] == "#":
            continue
        fields = line.split("\t", 8)
        if fields[2] == "gene" and 'gene_type "protein_coding"' in fields[8]:
            n_gene_total += 1
        elif fields[2] == "transcript" and 'transcript_type "protein_coding"' in fields[8]:
            n_tx_total += 1


def retained_frac_with_ci(gct, n_total, by_gene):
    """Proportion of the n_total-unit universe retained at each threshold, + 95% bootstrap
    CI over clusters. Universe = GCT rows plus (n_total - rows) never-expressed units."""
    df = pd.read_csv(gct, sep="\t", skiprows=2)
    mx = df[list(df.columns[2:])].max(axis=1).values
    ret = (mx[:, None] >= THRESHOLDS[None, :]).astype(float)  # units x thresholds
    cl = pd.factorize(df["Description"])[0] if by_gene else np.arange(len(df))
    n_clu = cl.max() + 1
    ret_c = np.zeros((n_clu, len(THRESHOLDS)))
    np.add.at(ret_c, cl, ret)
    size_c = np.bincount(cl, minlength=n_clu).astype(float)
    pad = n_total - len(df)  # singleton clusters of never-retained units
    ret_c = np.vstack([ret_c, np.zeros((pad, len(THRESHOLDS)))])
    size_c = np.concatenate([size_c, np.ones(pad)])
    n_all = len(size_c)
    counts = np.stack([np.bincount(rng.integers(0, n_all, n_all), minlength=n_all) for _ in range(N_BOOT)]).astype(float)
    boot = (counts @ ret_c) / (counts @ size_c)[:, None]
    point = ret_c.sum(axis=0) / n_total
    return point, np.percentile(boot, 2.5, axis=0), np.percentile(boot, 97.5, axis=0)


frac = {"gene": retained_frac_with_ci(IN_GENE_GCT, n_gene_total, by_gene=False),
        "transcript": retained_frac_with_ci(IN_TX_GCT, n_tx_total, by_gene=True)}
for k, (pt, lo, hi) in frac.items():
    print(f"{k} T=10: {pt[9]:.4f} [{lo[9]:.4f}, {hi[9]:.4f}]", file=sys.stderr)
dist = pd.read_csv(IN_DIST_CSV)

fig, (ax_ret, ax_dist) = plt.subplots(2, 1, figsize=(7.8, 16))  # stacked; panels 0.7x the old width
for name, src in [("gene", GENE_SRC), ("transcript", TX_SRC)]:
    sub = dist[(dist["source"] == src) & (dist["tpm_floor"] >= 1)]
    ax_dist.fill_between(sub["tpm_floor"], sub["ci_lo"], sub["ci_hi"], color=COLORS[name], alpha=0.25)
    ax_dist.plot(sub["tpm_floor"], sub["mean_dist"], marker="o", markersize=MARKER_SIZE, linewidth=LINE_WIDTH,
                 color=COLORS[name], label=name)
    point, lo, hi = frac[name]
    ax_ret.fill_between(THRESHOLDS, lo, hi, color=COLORS[name], alpha=0.25)
    ax_ret.plot(THRESHOLDS, point, marker="o", markersize=MARKER_SIZE, linewidth=LINE_WIDTH, color=COLORS[name], label=name)

ax_dist.set_ylabel("mean pairwise tissue distance (1 - r)", fontsize=AXIS_LABEL_FONTSIZE)
ax_ret.set_ylabel("proportion of all protein-coding\ngenes / transcripts retained", fontsize=AXIS_LABEL_FONTSIZE)
ax_ret.set_ylim(0, 1)
for ax in (ax_ret, ax_dist):
    ax.set_xlabel("TPM threshold", fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_xticks([1, 5, 10, 15, 20, 25, 30])
    ax.tick_params(labelsize=TICK_LABEL_FONTSIZE)
ax_ret.legend(fontsize=LEGEND_FONTSIZE, loc="upper right")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=600)
plt.close(fig)
print(f"gene_total={n_gene_total} tx_total={n_tx_total} -> {OUT_PNG}", file=sys.stderr)
