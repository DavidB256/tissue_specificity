"""Restrict a transcript-level GCT (already filtered to transcripts of
protein-coding GENES, e.g. GTEX_V11_TX_BASE_GCT) down further to transcripts
whose OWN GENCODE v47 transcript_type=="protein_coding" -- a strictly
narrower set, excluding e.g. nonsense_mediated_decay/retained_intron/
processed_transcript isoforms of an otherwise protein-coding gene. Base-ID
(version-stripped) match on transcript_id, same rationale as
filter_protein_coding_by_gene_id.py: GTEx pins different GENCODE ID versions
than v47."""
import sys
import re
import gzip
import pandas as pd

IN_GCT = sys.argv[1]
GENCODE_GTF = sys.argv[2]
OUT_GCT = sys.argv[3]

df = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(df.columns[2:])

pc_tx_ids = set()
with open(GENCODE_GTF) as f:
    for line in f:
        if line[0] == "#":
            continue
        fields = line.split("\t", 8)
        if fields[2] != "transcript":
            continue
        if 'transcript_type "protein_coding"' not in fields[8]:
            continue
        m = re.search(r'transcript_id "([^"]+)"', fields[8])
        pc_tx_ids.add(m.group(1).split(".")[0])

tx_base = df["Name"].str.split(".").str[0]
out = df[tx_base.isin(pc_tx_ids)]

print(f"GENCODE protein_coding transcript_ids: {len(pc_tx_ids)}", file=sys.stderr)
print(f"transcripts: {len(df)} -> {len(out)} protein_coding (by transcript_type)", file=sys.stderr)

with gzip.open(OUT_GCT, "wt") as f:
    f.write("#1.2\n")
    f.write(f"{len(out)}\t{len(tissues)}\n")
    out.to_csv(f, sep="\t", index=False)
