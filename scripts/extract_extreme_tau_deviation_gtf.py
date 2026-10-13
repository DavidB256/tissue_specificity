"""Filtered copy of the GENCODE GTF for IGV inspection: keep only transcripts
with arctanh(tx_tau_weighted) - arctanh(gene_tau_weighted) >= 3 AND
isoform_fraction (tx_total/gene_total TPM, summed across all tissues)
>= 10%, plus their exon/CDS records and their parent gene record. tau
scaled by 0.9999 before arctanh (arctanh(1) undefined). Sibling transcripts
of the same gene that don't meet the criteria are NOT included."""
import sys
import re
import pandas as pd
import numpy as np

IN_GTF = sys.argv[1]
IN_GCT = sys.argv[2]
IN_TX_TAU_CSV = sys.argv[3]
IN_INFO_CONTENT_CSV = sys.argv[4]
OUT_GTF = sys.argv[5]

SCALE = 0.9999
DEVIATION_THRESHOLD = 3.0
ISOFORM_FRACTION_THRESHOLD = 0.10


def arctanh_tau(tau):
    return np.arctanh(SCALE * tau)


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

tx_tau = pd.read_csv(IN_TX_TAU_CSV).rename(columns={
    "gene": "transcript_id", "tau_free_weighted": "tx_tau_weighted",
})[["transcript_id", "tx_tau_weighted"]]
info = pd.read_csv(IN_INFO_CONTENT_CSV)[["gene_id", "gene_tau_weighted"]]

df = tx.merge(tx_tau, on="transcript_id", how="inner")
df = df.merge(info, on="gene_id", how="inner")
df = df.dropna(subset=["tx_tau_weighted", "gene_tau_weighted", "isoform_fraction"])
df["deviation_arctanh"] = arctanh_tau(df["tx_tau_weighted"]) - arctanh_tau(df["gene_tau_weighted"])

selected = df[(df["deviation_arctanh"] >= DEVIATION_THRESHOLD) &
              (df["isoform_fraction"] >= ISOFORM_FRACTION_THRESHOLD)]
selected_tx = set(selected["transcript_id"])
selected_genes = set(selected["gene_id"])

n_written = {"gene": 0, "transcript": 0, "exon": 0, "CDS": 0}
with open(IN_GTF) as fin, open(OUT_GTF, "w") as fout:
    for line in fin:
        if line[0] == "#":
            fout.write(line)
            continue
        fields = line.split("\t", 8)
        feature = fields[2]
        attrs = fields[8]

        if feature == "gene":
            gene_id = re.search(r'gene_id "([^"]+)"', attrs).group(1)
            if gene_id in selected_genes:
                fout.write(line)
                n_written["gene"] += 1
        elif feature in ("transcript", "exon", "CDS"):
            transcript_id = re.search(r'transcript_id "([^"]+)"', attrs).group(1)
            if transcript_id in selected_tx:
                fout.write(line)
                n_written[feature] += 1

selected.to_csv(OUT_GTF + ".selected_tx.csv", index=False)

print(f"n_selected_tx={len(selected_tx)} n_selected_genes={len(selected_genes)} "
      f"lines_written={n_written} -> {OUT_GTF}", file=sys.stderr)
