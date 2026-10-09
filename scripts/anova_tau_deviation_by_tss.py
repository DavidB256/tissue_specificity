"""Test hypothesis: TSS choice drives variation in tx tau_weighted among a
gene's isoforms. For each multi-tx gene, group its quantified transcripts
(tx_total>0) into categories by exact TSS coordinate (transcript_tss.csv).
One-way ANOVA on (tx_tau_weighted - gene_tau_weighted) across TSS groups.
Only genes with >=2 distinct TSS groups and N>k (nonzero error df, i.e. at
least one TSS group shared by >=2 transcripts) are testable."""
import sys
import numpy as np
import pandas as pd
from scipy.stats import f_oneway

IN_GCT = sys.argv[1]
IN_TX_TAU_CSV = sys.argv[2]
IN_INFO_CONTENT_CSV = sys.argv[3]
IN_TSS_CSV = sys.argv[4]
OUT_CSV = sys.argv[5]
OUT_PNG = sys.argv[6]

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

rows = []
for gene_id, g in df.groupby("gene_id"):
    if len(g) < 3:
        continue
    groups = [grp["deviation"].values for _, grp in g.groupby("tss")]
    k = len(groups)
    n = len(g)
    if k < 2 or n <= k:
        continue
    f_stat, p_value = f_oneway(*groups)
    if np.isnan(p_value):
        continue
    rows.append({
        "gene_id": gene_id,
        "gene_tau_weighted": g["gene_tau_weighted"].iloc[0],
        "n_tx": n,
        "n_tss_groups": k,
        "f_stat": f_stat,
        "p_value": p_value,
    })

out = pd.DataFrame(rows).sort_values("p_value")
out.to_csv(OUT_CSV, index=False)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

neg_log10_p = -np.log10(out["p_value"].clip(lower=1e-300))

fig, ax = plt.subplots(figsize=(8, 7))
ax.scatter(out["gene_tau_weighted"], neg_log10_p, s=5, alpha=0.3)
ax.axhline(-np.log10(0.05), color="gray", linestyle="--", linewidth=1, label="p=0.05")
ax.set_xlabel("gene tau_weighted")
ax.set_ylabel("-log10(ANOVA p-value)")
ax.set_title(
    f"gtex_v11_tx_protein_coding: tau deviation ANOVA by TSS group\n"
    f"n_genes_tested={len(out)}, n_p<0.05={int((out['p_value'] < 0.05).sum())}"
)
ax.legend()
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)

print(f"n_genes_tested={len(out)} n_p<0.05={(out['p_value'] < 0.05).sum()} -> {OUT_CSV}, {OUT_PNG}",
      file=sys.stderr)
