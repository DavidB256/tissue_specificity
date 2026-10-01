"""Restrict a transcript-level GCT (Name=transcript_id) to protein-coding transcripts,
per gene_type=="protein_coding" in a transcript metadata CSV (transcript_id,gene_type
-- build_chess_transcript_meta.py). Generic (not GENCODE-specific like
filter_protein_coding.py): exact Name match against the meta table, no base-ID
stripping needed since both sides are the same CHESS transcript IDs."""
import sys
import gzip
import pandas as pd

IN_GCT = sys.argv[1]
META_CSV = sys.argv[2]
OUT_GCT = sys.argv[3]

df = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(df.columns[2:])

gene_type = pd.read_csv(META_CSV).set_index("transcript_id")["gene_type"]
is_pc = df["Name"].map(gene_type) == "protein_coding"
out = df[is_pc]

print(f"transcripts: {len(df)} -> {len(out)} protein_coding", file=sys.stderr)

with gzip.open(OUT_GCT, "wt") as f:
    f.write("#1.2\n")
    f.write(f"{len(out)}\t{len(tissues)}\n")
    out.to_csv(f, sep="\t", index=False)
