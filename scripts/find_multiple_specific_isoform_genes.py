""""multiple specific isoform" genes: genes with >=2 transcripts that each
satisfy (i) isoform_fraction (tx_total/gene_total TPM, summed across all
tissues) >= 1%, (ii) prob_tx_tau_gt_gene_weighted >= 0.95 (donor-level
bootstrap: tx tau_weighted exceeds gene tau_weighted in >=95% of
resamples)."""
import sys
import pandas as pd

IN_GCT = sys.argv[1]
IN_BOOTSTRAP_TX_TAU_CSV = sys.argv[2]
OUT_CSV = sys.argv[3]

ISOFORM_FRACTION_THRESHOLD = 0.01
PROB_THRESHOLD = 0.95

gct = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(gct.columns[2:])
tx_total = gct[tissues].sum(axis=1)
keep = tx_total > 0

tx = pd.DataFrame({
    "transcript_id": gct["Name"][keep],
    "gene_id": gct["Description"][keep],
    "tx_total": tx_total[keep],
})
tx["gene_total"] = tx.groupby("gene_id")["tx_total"].transform("sum")
tx["isoform_fraction"] = tx["tx_total"] / tx["gene_total"]

boot = pd.read_csv(IN_BOOTSTRAP_TX_TAU_CSV)[["transcript_id", "gene_id", "prob_tx_tau_gt_gene_weighted"]]

df = tx.merge(boot, on=["transcript_id", "gene_id"], how="inner")
df = df.dropna(subset=["isoform_fraction", "prob_tx_tau_gt_gene_weighted"])

qualifying = df[(df["isoform_fraction"] >= ISOFORM_FRACTION_THRESHOLD) &
                (df["prob_tx_tau_gt_gene_weighted"] >= PROB_THRESHOLD)]

n_qualifying = qualifying.groupby("gene_id").size().rename("n_qualifying_tx")
multi = n_qualifying[n_qualifying >= 2]

qualifying_tx_lists = (qualifying[qualifying["gene_id"].isin(multi.index)]
                        .groupby("gene_id")["transcript_id"].apply(lambda s: ",".join(sorted(s))))

out = pd.DataFrame({
    "gene_id": multi.index,
    "n_qualifying_tx": multi.values,
    "qualifying_transcript_ids": qualifying_tx_lists.reindex(multi.index).values,
}).sort_values("n_qualifying_tx", ascending=False)
out.to_csv(OUT_CSV, index=False)

print(f"n_genes={len(out)} -> {OUT_CSV}", file=sys.stderr)
