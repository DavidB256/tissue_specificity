"""transcript_id,gene_id,gene_name table for chess_tx_protein_coding -- transcript_id/
gene_id read straight off the (already protein_coding-filtered) GCT's Name/Description
columns, gene_name joined from IN_FULL_META_CSV (parse_chess_gtf_meta.py) by
transcript_id. Parallel in shape to build_gtex_v11_tx_meta.py's output."""
import sys
import pandas as pd

IN_GCT = sys.argv[1]
IN_FULL_META_CSV = sys.argv[2]
OUT_CSV = sys.argv[3]

gct = pd.read_csv(IN_GCT, sep="\t", skiprows=2, usecols=["Name", "Description"])
df = pd.DataFrame({"transcript_id": gct["Name"], "gene_id": gct["Description"]})

full_meta = pd.read_csv(IN_FULL_META_CSV)[["transcript_id", "gene_name"]]
df = df.merge(full_meta, on="transcript_id", how="left")

df.to_csv(OUT_CSV, index=False)

n_unmapped = df["gene_name"].isna().sum()
print(f"transcripts={len(df)} unmapped_gene_name={n_unmapped} -> {OUT_CSV}", file=sys.stderr)
