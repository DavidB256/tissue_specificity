"""Restrict a transcript-level GCT keyed by Ensembl gene_id in Description (e.g.
gtex_v11_tx, RSEM/GTEx's own gene_id, not necessarily GENCODE v47's exact ID version)
to protein-coding transcripts, per GENCODE v47 gene-level gene_type=="protein_coding".
Base-ID (version-stripped) match, same rationale as filter_protein_coding.py: GTEx
pins different GENCODE ID versions than v47, so an exact versioned match would miss
most genes."""
import sys
import re
import gzip
import pandas as pd

IN_GCT = sys.argv[1]
GENCODE_GTF = sys.argv[2]
OUT_GCT = sys.argv[3]

df = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(df.columns[2:])

pc_gene_ids = set()
with open(GENCODE_GTF) as f:
    for line in f:
        if line[0] == "#":
            continue
        fields = line.split("\t", 8)
        if fields[2] != "gene":
            continue
        if 'gene_type "protein_coding"' not in fields[8]:
            continue
        m = re.search(r'gene_id "([^"]+)"', fields[8])
        pc_gene_ids.add(m.group(1).split(".")[0])

gene_base = df["Description"].str.split(".").str[0]
out = df[gene_base.isin(pc_gene_ids)]

print(f"GENCODE protein_coding gene_ids: {len(pc_gene_ids)}", file=sys.stderr)
print(f"transcripts: {len(df)} -> {len(out)} protein_coding", file=sys.stderr)

with gzip.open(OUT_GCT, "wt") as f:
    f.write("#1.2\n")
    f.write(f"{len(out)}\t{len(tissues)}\n")
    out.to_csv(f, sep="\t", index=False)
