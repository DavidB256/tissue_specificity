"""2-panel violin fig: x = number of distinct TSSs among a gene's
"significant" transcripts, y = variance of tx_tau_weighted among those
transcripts. Panel 1: significant = isoform_fraction >= 10%. Panel 2: >= 1%.
Genes need >=2 significant tx (variance undefined otherwise). IQR drawn per
violin."""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IN_GCT = sys.argv[1]
IN_TX_TAU_CSV = sys.argv[2]
IN_TSS_CSV = sys.argv[3]
OUT_PNG = sys.argv[4]

THRESHOLDS = [0.10, 0.01]

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

tx_tau = pd.read_csv(IN_TX_TAU_CSV).rename(columns={
    "gene": "transcript_id", "tau_free_weighted": "tx_tau_weighted",
})[["transcript_id", "tx_tau_weighted"]]
tss = pd.read_csv(IN_TSS_CSV)[["transcript_id", "tss"]]

df = isoform_fraction_df.merge(tx_tau, on="transcript_id", how="inner")
df = df.merge(tss, on="transcript_id", how="inner")
df = df.dropna(subset=["isoform_fraction", "tx_tau_weighted", "tss"])


def gene_stats(threshold):
    sig = df[df["isoform_fraction"] >= threshold]
    stats = sig.groupby("gene_id").agg(
        n_sig_tx=("tx_tau_weighted", "count"),
        n_distinct_tss=("tss", "nunique"),
        tau_var=("tx_tau_weighted", "var"),
    )
    return stats[stats["n_sig_tx"] >= 2]


fig, axes = plt.subplots(1, 2, figsize=(15, 7))
for ax, threshold in zip(axes, THRESHOLDS):
    stats = gene_stats(threshold)
    tss_counts = sorted(stats["n_distinct_tss"].unique())
    groups = [stats.loc[stats["n_distinct_tss"] == c, "tau_var"].dropna().values for c in tss_counts]
    ns = [len(g) for g in groups]

    parts = ax.violinplot(groups, positions=tss_counts, showmedians=True, widths=0.7,
                           quantiles=[[0.25, 0.75]] * len(groups))
    for pc in parts["bodies"]:
        pc.set_facecolor("steelblue")
        pc.set_alpha(0.5)

    ax.set_xticks(tss_counts)
    ax.set_xticklabels([f"{c}\nn={n}" for c, n in zip(tss_counts, ns)], fontsize=8)
    ax.set_xlabel("n distinct TSS among significant tx")
    ax.set_ylabel("variance of tx tau_weighted")
    ax.set_title(f"significant tx = isoform_fraction >= {threshold:.0%}\nn_genes={sum(ns)}")

    print(f"threshold={threshold} n_genes={sum(ns)} tss_counts={tss_counts} ns={ns}", file=sys.stderr)

fig.suptitle("gtex_v11_tx_protein_coding: tau_weighted variance vs. n distinct TSS, by significance threshold")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)
