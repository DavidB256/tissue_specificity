"""Poster fig, 2 panels vs. TPM threshold T (x: 1..30; retained = median TPM >= T in at
least one tissue): left, mean pairwise tissue distance (1 - r) with 95% bootstrap band
(precomputed by plot_mean_tissue_dist_vs_tpm_floor.py), gene vs. transcript; right,
proportion of all protein-coding genes / protein-coding transcripts in GENCODE v47 that
are retained. Colors: matplotlib defaults (gene blue, transcript orange). dpi=600.
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


def retained_counts(gct):
    df = pd.read_csv(gct, sep="\t", skiprows=2)
    mx = df[list(df.columns[2:])].max(axis=1).values
    return np.array([(mx >= t).sum() for t in THRESHOLDS])


frac = {"gene": retained_counts(IN_GENE_GCT) / n_gene_total,
        "transcript": retained_counts(IN_TX_GCT) / n_tx_total}
dist = pd.read_csv(IN_DIST_CSV)

fig, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(19, 8))
for name, src in [("gene", GENE_SRC), ("transcript", TX_SRC)]:
    sub = dist[(dist["source"] == src) & (dist["tpm_floor"] >= 1)]
    ax_l.fill_between(sub["tpm_floor"], sub["ci_lo"], sub["ci_hi"], color=COLORS[name], alpha=0.25)
    ax_l.plot(sub["tpm_floor"], sub["mean_dist"], marker="o", color=COLORS[name], label=name)
    ax_r.plot(THRESHOLDS, frac[name], marker="o", color=COLORS[name], label=name)

ax_l.set_ylabel("mean pairwise tissue distance (1 - r)", fontsize=AXIS_LABEL_FONTSIZE)
ax_r.set_ylabel("proportion of all protein-coding\ngenes / transcripts retained", fontsize=AXIS_LABEL_FONTSIZE)
ax_r.set_ylim(0, 1)
for ax in (ax_l, ax_r):
    ax.set_xlabel("TPM threshold", fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_xticks([1, 5, 10, 15, 20, 25, 30])
    ax.tick_params(labelsize=TICK_LABEL_FONTSIZE)
ax_l.legend(fontsize=LEGEND_FONTSIZE, loc="upper left")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=600)
plt.close(fig)
print(f"gene_total={n_gene_total} tx_total={n_tx_total} -> {OUT_PNG}", file=sys.stderr)
