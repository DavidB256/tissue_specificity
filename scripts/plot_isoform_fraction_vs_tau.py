"""2-row, 3-col scatter, one point per transcript (multi-tx genes only, n_tx>=2),
x = isoform_fraction (tx_total / gene_total TPM, from the base
gtex_v11_tx_protein_coding GCT, tx_total>0 only).

Row 1 (tau_weighted, tau_free_weighted):
  (1) y = tx tau_weighted
  (2) y = tx tau_weighted - gene tau_weighted
  (3) y = prob_tx_tau_gt_gene_weighted (bootstrap, proportion of resamples in
      which the tx's tau_weighted exceeds its gene's tau_weighted)
Row 2 (tau_unweighted): same 3 quantities, unweighted tau instead.

Points are colored by whether the two-sided bootstrap p-value for "tx tau >
gene tau" (p = 2*min(prob, 1-prob)) is <0.05 -- row 1 uses
prob_tx_tau_gt_gene_weighted, row 2 uses prob_tx_tau_gt_gene_unweighted."""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from gls_tau import weighted_tau

IN_GCT = sys.argv[1]
IN_TX_TAU_CSV = sys.argv[2]
IN_INFO_CONTENT_CSV = sys.argv[3]
IN_BOOTSTRAP_TX_TAU_CSV = sys.argv[4]
OUT_PNG = sys.argv[5]

gct = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(gct.columns[2:])
tx_total = gct[tissues].sum(axis=1)
keep = tx_total > 0

isoform_fraction_df = pd.DataFrame({
    "transcript_id": gct["Name"][keep],
    "gene_id": gct["Description"][keep],
    "tx_total": tx_total[keep],
})
isoform_fraction_df["gene_total"] = isoform_fraction_df.groupby("gene_id")["tx_total"].transform("sum")
isoform_fraction_df["isoform_fraction"] = isoform_fraction_df["tx_total"] / isoform_fraction_df["gene_total"]

# gene tau_unweighted: not precomputed anywhere for this exact (PC, tx_total>0) universe
# -- recompute directly from the same GCT so it matches the isoform_fraction universe.
gene_mat = gct.loc[keep, tissues].groupby(gct.loc[keep, "Description"]).sum()
gene_tau_unweighted = pd.Series(
    weighted_tau(gene_mat.to_numpy(), np.ones(len(tissues))),
    index=gene_mat.index, name="gene_tau_unweighted",
)

tx_tau = pd.read_csv(IN_TX_TAU_CSV).rename(columns={
    "gene": "transcript_id", "tau_unweighted": "tx_tau_unweighted", "tau_free_weighted": "tx_tau_weighted",
})
info = pd.read_csv(IN_INFO_CONTENT_CSV)[["gene_id", "n_tx", "gene_tau_weighted"]]
boot = pd.read_csv(IN_BOOTSTRAP_TX_TAU_CSV)[
    ["transcript_id", "prob_tx_tau_gt_gene_weighted", "prob_tx_tau_gt_gene_unweighted"]
]

df = isoform_fraction_df.merge(tx_tau, on="transcript_id", how="inner")
df = df.merge(info, on="gene_id", how="inner")
df = df.merge(gene_tau_unweighted, left_on="gene_id", right_index=True, how="inner")
df = df.merge(boot, on="transcript_id", how="inner")
df = df[df["n_tx"] >= 2]

df["tau_diff_weighted"] = df["tx_tau_weighted"] - df["gene_tau_weighted"]
df["tau_diff_unweighted"] = df["tx_tau_unweighted"] - df["gene_tau_unweighted"]
df["p_weighted"] = 2 * np.minimum(df["prob_tx_tau_gt_gene_weighted"], 1 - df["prob_tx_tau_gt_gene_weighted"])
df["p_unweighted"] = 2 * np.minimum(df["prob_tx_tau_gt_gene_unweighted"], 1 - df["prob_tx_tau_gt_gene_unweighted"])

df = df.dropna(subset=[
    "isoform_fraction", "tx_tau_weighted", "tau_diff_weighted", "prob_tx_tau_gt_gene_weighted",
    "tx_tau_unweighted", "tau_diff_unweighted", "prob_tx_tau_gt_gene_unweighted",
])

SIG_COLOR, NONSIG_COLOR = "crimson", "steelblue"


def scatter_by_sig(ax, x, y, sig):
    ax.scatter(x[~sig], y[~sig], s=8, alpha=0.1, color=NONSIG_COLOR, label="p>=0.05")
    ax.scatter(x[sig], y[sig], s=8, alpha=0.1, color=SIG_COLOR, label="p<0.05")


rows = [
    ("weighted", "p_weighted", "tx_tau_weighted", "tau_diff_weighted", "prob_tx_tau_gt_gene_weighted",
     "tx tau_weighted", "tx tau_weighted - gene tau_weighted", "P(tx tau_weighted > gene tau_weighted)"),
    ("unweighted", "p_unweighted", "tx_tau_unweighted", "tau_diff_unweighted", "prob_tx_tau_gt_gene_unweighted",
     "tx tau_unweighted", "tx tau_unweighted - gene tau_unweighted", "P(tx tau_unweighted > gene tau_unweighted)"),
]

fig, axes = plt.subplots(2, 3, figsize=(18, 12))
for r, (tag, pcol, y1, y2, y3, ylab1, ylab2, ylab3) in enumerate(rows):
    sig = df[pcol] < 0.05
    for c, (ycol, ylab) in enumerate([(y1, ylab1), (y2, ylab2), (y3, ylab3)]):
        ax = axes[r, c]
        scatter_by_sig(ax, df["isoform_fraction"], df[ycol], sig)
        ax.set_xlabel("isoform fraction")
        ax.set_ylabel(ylab)
        if c == 1:
            ax.axhline(0, color="black", linestyle="--", linewidth=1)
    axes[r, 0].legend(loc="best", fontsize=8, markerscale=3)

fig.suptitle(f"gtex_v11_tx_protein_coding: isoform fraction vs. tau, multi-tx genes only (n={len(df)})")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)

print(f"n={len(df)} n_sig_weighted={int((df['p_weighted']<0.05).sum())} "
      f"n_sig_unweighted={int((df['p_unweighted']<0.05).sum())} -> {OUT_PNG}", file=sys.stderr)
