"""Restrict the Baron GCT (genes keyed by symbol, no Ensembl IDs -- see
build_baron_gct.py) to protein-coding genes only, per GENCODE gene_type==
"protein_coding" gene-level rows, matched on gene_name (symbol). A handful of
symbols (8 in GENCODE v47) are used for both a protein_coding gene and a
non-coding gene at different loci; treated as protein_coding if any matching
GENCODE gene row is protein_coding (permissive, mirrors filter_protein_coding.py's
GTEx behavior of trusting GENCODE's own gene_type field)."""
import sys
import re
import gzip
import pandas as pd

IN_GCT = sys.argv[1]
GENCODE_GTF = sys.argv[2]
OUT_GCT = sys.argv[3]

pc_symbols = set()
with open(GENCODE_GTF) as f:
    for line in f:
        if line[0] == "#":
            continue
        fields = line.split("\t", 8)
        if fields[2] != "gene":
            continue
        if 'gene_type "protein_coding"' not in fields[8]:
            continue
        gene_name = re.search(r'gene_name "([^"]+)"', fields[8]).group(1)
        pc_symbols.add(gene_name)

df = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
out = df[df["Name"].isin(pc_symbols)]
cell_types = list(df.columns[2:])

print(f"GENCODE protein_coding symbols: {len(pc_symbols)}", file=sys.stderr)
print(f"genes: {len(df)} -> {len(out)} protein_coding "
      f"({len(df) - len(out)} dropped: not a GENCODE protein_coding symbol)", file=sys.stderr)

with gzip.open(OUT_GCT, "wt") as f:
    f.write("#1.2\n")
    f.write(f"{len(out)}\t{len(cell_types)}\n")
    out.to_csv(f, sep="\t", index=False)
