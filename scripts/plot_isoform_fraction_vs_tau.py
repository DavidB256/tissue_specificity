"""3-panel scatter, one point per transcript (multi-tx genes only, n_tx>=2),
x = isoform_fraction (tx_total / gene_total TPM, from the base
gtex_v11_tx_protein_coding GCT, tx_total>0 only):
  (1) y = tx tau_free_weighted
  (2) y = tx tau_free_weighted - gene_tau_weighted
  (3) y = prob_tx_tau_gt_gene_weighted (bootstrap, proportion of resamples in
      which the tx's tau_weighted exceeds its gene's tau_weighted)"""
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

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

tx_tau = pd.read_csv(IN_TX_TAU_CSV).rename(columns={"gene": "transcript_id", "tau_free_weighted": "tx_tau_weighted"})
info = pd.read_csv(IN_INFO_CONTENT_CSV)[["gene_id", "n_tx", "gene_tau_weighted"]]
boot = pd.read_csv(IN_BOOTSTRAP_TX_TAU_CSV)[["transcript_id", "prob_tx_tau_gt_gene_weighted"]]

df = isoform_fraction_df.merge(tx_tau[["transcript_id", "tx_tau_weighted"]], on="transcript_id", how="inner")
df = df.merge(info, on="gene_id", how="inner")
df = df.merge(boot, on="transcript_id", how="inner")
df = df[df["n_tx"] >= 2]
df["tau_diff"] = df["tx_tau_weighted"] - df["gene_tau_weighted"]
df = df.dropna(subset=["isoform_fraction", "tx_tau_weighted", "tau_diff", "prob_tx_tau_gt_gene_weighted"])

fig, axes = plt.subplots(1, 3, figsize=(18, 6))

axes[0].scatter(df["isoform_fraction"], df["tx_tau_weighted"], s=8, alpha=0.1, color="steelblue")
axes[0].set_xlabel("isoform fraction")
axes[0].set_ylabel("tx tau_weighted")

axes[1].scatter(df["isoform_fraction"], df["tau_diff"], s=8, alpha=0.1, color="steelblue")
axes[1].axhline(0, color="black", linestyle="--", linewidth=1)
axes[1].set_xlabel("isoform fraction")
axes[1].set_ylabel("tx tau_weighted - gene tau_weighted")

axes[2].scatter(df["isoform_fraction"], df["prob_tx_tau_gt_gene_weighted"], s=8, alpha=0.1, color="steelblue")
axes[2].set_xlabel("isoform fraction")
axes[2].set_ylabel("P(tx tau_weighted > gene tau_weighted)")

fig.suptitle(f"gtex_v11_tx_protein_coding: isoform fraction vs. tau_weighted, multi-tx genes only (n={len(df)})")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)

print(f"n={len(df)} -> {OUT_PNG}", file=sys.stderr)
