""""specific_isoform" genes (gtex_v11_tx_protein_coding, multi-tx genes only,
n_tx>=2): genes where (i) a subset of transcripts has
prob_tx_tau_gt_gene_weighted >= 0.90 (bootstrap: tx tau_weighted exceeds gene
tau_weighted in >=90% of resamples) and (ii) that subset together composes
>=90% of the gene's total expression (sum of isoform_fraction = tx_total /
gene_total, from the base GCT, tx_total>0 only)."""
import sys
import numpy as np
import pandas as pd
from gls_tau import weighted_tau

IN_GCT = sys.argv[1]
IN_BOOTSTRAP_TX_TAU_CSV = sys.argv[2]
IN_INFO_CONTENT_CSV = sys.argv[3]
OUT_CSV = sys.argv[4]

PROB_THRESHOLD = 0.90
FRACTION_THRESHOLD = 0.90

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

gene_mat = gct.loc[keep, tissues].groupby(gct.loc[keep, "Description"]).sum()
gene_tau_unweighted = pd.Series(
    weighted_tau(gene_mat.to_numpy(), np.ones(len(tissues))),
    index=gene_mat.index, name="gene_tau_unweighted",
)

boot = pd.read_csv(IN_BOOTSTRAP_TX_TAU_CSV)[["transcript_id", "gene_id", "prob_tx_tau_gt_gene_weighted"]]
info = pd.read_csv(IN_INFO_CONTENT_CSV)[["gene_id", "n_tx", "gene_tau_weighted"]]

df = boot.merge(isoform_fraction_df[["transcript_id", "isoform_fraction"]], on="transcript_id", how="inner")
df = df.merge(info, on="gene_id", how="inner")
df = df[df["n_tx"] >= 2]

df["high_tau"] = df["prob_tx_tau_gt_gene_weighted"] >= PROB_THRESHOLD
frac_high_tau = (
    df[df["high_tau"]].groupby("gene_id")["isoform_fraction"].sum()
    .reindex(df["gene_id"].unique(), fill_value=0.0)
    .rename("isoform_fraction_with_high_tau")
)

genes = df[["gene_id", "n_tx", "gene_tau_weighted"]].drop_duplicates("gene_id").set_index("gene_id")
genes = genes.join(frac_high_tau).join(gene_tau_unweighted)

specific_isoform = genes[genes["isoform_fraction_with_high_tau"] >= FRACTION_THRESHOLD]
out = specific_isoform.reset_index()[
    ["gene_id", "isoform_fraction_with_high_tau", "gene_tau_unweighted", "gene_tau_weighted"]
].sort_values("isoform_fraction_with_high_tau", ascending=False)
out.to_csv(OUT_CSV, index=False)

print(f"multi-tx genes considered: {len(genes)}", file=sys.stderr)
print(f"specific_isoform genes (frac_high_tau>={FRACTION_THRESHOLD}): {len(out)} -> {OUT_CSV}", file=sys.stderr)
