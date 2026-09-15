"""Restrict a GCT to protein-coding genes only, per GENCODE gene_type=="protein_coding"
(gene-level rows of the GENCODE annotation GTF). Auto-detects the GCT's gene key from
the "Name" column: Ensembl IDs (GTEx, e.g. "ENSG00000141510.16") match on base ID
(version suffix stripped) -- exact versioned-ID match against GENCODE v47 covers <4%
of GTEx's genes (GTEx v11 pins different GENCODE ID versions than v47), while base-ID
match recovers ~96% of GENCODE's protein_coding set. Gene symbols (Baron, no Ensembl
IDs -- see build_baron_gct.py) match on gene_name directly; a handful of symbols (8 in
GENCODE v47) are used for both a protein_coding gene and a non-coding gene at different
loci, treated as protein_coding if any matching GENCODE gene row is protein_coding
(permissive, trusts GENCODE's own gene_type field either way). Rows are kept in their
original order with their original IDs untouched -- only whole rows are dropped."""
import sys
import re
import gzip
import pandas as pd

IN_GCT = sys.argv[1]
GENCODE_GTF = sys.argv[2]
OUT_GCT = sys.argv[3]

df = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
is_ensembl = df["Name"].str.startswith("ENSG").mean() > 0.5
key = "gene_id" if is_ensembl else "gene_name"

pc_keys = set()
with open(GENCODE_GTF) as f:
    for line in f:
        if line[0] == "#":
            continue
        fields = line.split("\t", 8)
        if fields[2] != "gene":
            continue
        if 'gene_type "protein_coding"' not in fields[8]:
            continue
        value = re.search(rf'{key} "([^"]+)"', fields[8]).group(1)
        pc_keys.add(value.split(".")[0] if is_ensembl else value)

gene_key = df["Name"].str.split(".").str[0] if is_ensembl else df["Name"]
out = df[gene_key.isin(pc_keys)]
tissues = list(df.columns[2:])

print(f"key={key} GENCODE protein_coding {key}s: {len(pc_keys)}", file=sys.stderr)
print(f"genes: {len(df)} -> {len(out)} protein_coding "
      f"({len(pc_keys) - out.shape[0]} GENCODE protein_coding {key}s not found in GCT)", file=sys.stderr)

with gzip.open(OUT_GCT, "wt") as f:
    f.write("#1.2\n")
    f.write(f"{len(out)}\t{len(tissues)}\n")
    out.to_csv(f, sep="\t", index=False)
