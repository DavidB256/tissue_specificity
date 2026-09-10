"""Restrict a GCT to protein-coding genes only, per GENCODE gene_type=="protein_coding"
(gene-level rows of the GENCODE annotation GTF). Matched on base Ensembl gene ID (version
suffix stripped): exact versioned-ID match against GENCODE v47 covers <4% of GTEx's genes
(GTEx v11 pins different GENCODE ID versions than v47), while base-ID match recovers
~96% of GENCODE's protein_coding set. Rows are kept in their original order with their
original (versioned) IDs untouched -- only whole rows are dropped."""
import sys
import re
import gzip
import pandas as pd

IN_GCT = sys.argv[1]
GENCODE_GTF = sys.argv[2]
OUT_GCT = sys.argv[3]

pc_base_ids = set()
with open(GENCODE_GTF) as f:
    for line in f:
        if line[0] == "#":
            continue
        fields = line.split("\t", 8)
        if fields[2] != "gene":
            continue
        if 'gene_type "protein_coding"' not in fields[8]:
            continue
        gene_id = re.search(r'gene_id "([^"]+)"', fields[8]).group(1)
        pc_base_ids.add(gene_id.split(".")[0])

df = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
base_id = df["Name"].str.split(".").str[0]
out = df[base_id.isin(pc_base_ids)]
tissues = list(df.columns[2:])

print(f"GENCODE protein_coding genes: {len(pc_base_ids)}", file=sys.stderr)
print(f"genes: {len(df)} -> {len(out)} protein_coding "
      f"({len(pc_base_ids) - out.shape[0]} GENCODE protein_coding genes not found in GCT)", file=sys.stderr)

with gzip.open(OUT_GCT, "wt") as f:
    f.write("#1.2\n")
    f.write(f"{len(out)}\t{len(tissues)}\n")
    out.to_csv(f, sep="\t", index=False)
