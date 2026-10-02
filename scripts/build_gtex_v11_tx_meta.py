"""transcript_id,gene_id,gene_name table for gtex_v11_tx -- transcript_id/gene_id read
straight off the (already protein_coding-filtered) GCT's Name/Description columns,
gene_name joined from gencode_v47_gene_id_name_map.csv on base (version-stripped)
gene_id, same rationale as filter_protein_coding_by_gene_id.py. Only these 3 columns
are needed: parallel in shape to chess_transcript_meta.csv but without CHESS-specific
gene_type/length_bp, since this source's protein_coding filter and TPM units don't
need either."""
import sys
import pandas as pd

IN_GCT = sys.argv[1]
IN_NAME_MAP_CSV = sys.argv[2]
OUT_CSV = sys.argv[3]

gct = pd.read_csv(IN_GCT, sep="\t", skiprows=2, usecols=["Name", "Description"])
df = pd.DataFrame({"transcript_id": gct["Name"], "gene_id": gct["Description"]})
df["gene_base"] = df["gene_id"].str.split(".").str[0]

name_map = pd.read_csv(IN_NAME_MAP_CSV)
name_map["gene_base"] = name_map["gene_id"].str.split(".").str[0]
df = df.merge(name_map[["gene_base", "gene_name"]], on="gene_base", how="left").drop(columns="gene_base")

df.to_csv(OUT_CSV, index=False)

n_unmapped = df["gene_name"].isna().sum()
print(f"transcripts={len(df)} unmapped_gene_name={n_unmapped} -> {OUT_CSV}", file=sys.stderr)
