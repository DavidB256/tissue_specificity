"""Pooled one-way-random-effects ICC: how much of the variance in
(tx_tau_weighted - gene_tau_weighted) sits between TSS-groups vs. within a
TSS-group, pooled across all multi-tx genes. Group = (gene_id, TSS-cluster),
TSSs within 50bp merged (see plot_tau_variance_by_tss_count_grouped50bp.py).
All quantified transcripts (tx_total>0) used, not just "significant" ones."""
import sys
import numpy as np
import pandas as pd

IN_GCT = sys.argv[1]
IN_TX_TAU_CSV = sys.argv[2]
IN_INFO_CONTENT_CSV = sys.argv[3]
IN_TSS_CSV = sys.argv[4]
OUT_CSV = sys.argv[5]

MAX_GAP = 50

gct = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(gct.columns[2:])
tx_total = gct[tissues].sum(axis=1)
keep = tx_total > 0

tx = pd.DataFrame({"transcript_id": gct["Name"][keep], "gene_id": gct["Description"][keep]})

tx_tau = pd.read_csv(IN_TX_TAU_CSV).rename(columns={
    "gene": "transcript_id", "tau_free_weighted": "tx_tau_weighted",
})[["transcript_id", "tx_tau_weighted"]]
info = pd.read_csv(IN_INFO_CONTENT_CSV)[["gene_id", "gene_tau_weighted"]]
tss = pd.read_csv(IN_TSS_CSV)[["transcript_id", "tss"]]

df = tx.merge(tx_tau, on="transcript_id", how="inner")
df = df.merge(info, on="gene_id", how="inner")
df = df.merge(tss, on="transcript_id", how="inner")
df = df.dropna(subset=["tx_tau_weighted", "gene_tau_weighted", "tss"])
df["deviation"] = df["tx_tau_weighted"] - df["gene_tau_weighted"]

n_tx_per_gene = df.groupby("gene_id")["transcript_id"].transform("count")
df = df[n_tx_per_gene >= 2].copy()


def assign_tss_clusters(tss_values):
    order = np.argsort(tss_values)
    sorted_vals = tss_values.values[order]
    cluster = np.zeros(len(sorted_vals), dtype=int)
    c = 0
    for i in range(1, len(sorted_vals)):
        if sorted_vals[i] - sorted_vals[i - 1] > MAX_GAP:
            c += 1
        cluster[i] = c
    out = np.empty(len(sorted_vals), dtype=int)
    out[order] = cluster
    return out


df["tss_cluster"] = df.groupby("gene_id")["tss"].transform(assign_tss_clusters)
df["group"] = list(zip(df["gene_id"], df["tss_cluster"]))

groups = df.groupby("group")["deviation"]
group_n = groups.size()
group_mean = groups.mean()

N = len(df)
k = len(group_n)
grand_mean = df["deviation"].mean()

ss_between = (group_n * (group_mean - grand_mean) ** 2).sum()
ms_between = ss_between / (k - 1)

ss_within = groups.apply(lambda s: ((s - s.mean()) ** 2).sum()).sum()
df_within = N - k
ms_within = ss_within / df_within if df_within > 0 else np.nan

n0 = (N - (group_n ** 2).sum() / N) / (k - 1)

icc = (ms_between - ms_within) / (ms_between + (n0 - 1) * ms_within)
icc_clipped = max(icc, 0.0)

out = pd.DataFrame([{
    "n_genes": df["gene_id"].nunique(),
    "n_tx": N,
    "n_groups": k,
    "n0": n0,
    "ms_between": ms_between,
    "ms_within": ms_within,
    "icc": icc,
    "icc_clipped": icc_clipped,
}])
out.to_csv(OUT_CSV, index=False)

print(f"n_genes={out['n_genes'].iloc[0]} n_tx={N} n_groups={k} n0={n0:.3f} "
      f"ms_between={ms_between:.6f} ms_within={ms_within:.6f} icc={icc:.4f} "
      f"icc_clipped={icc_clipped:.4f} -> {OUT_CSV}", file=sys.stderr)
