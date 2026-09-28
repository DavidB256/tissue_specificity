"""Map Ensembl gene ID (versioned, as listed in the GTF's gene_id attribute) to gene
name, from gene-level rows of a GENCODE annotation GTF. One row per gene_id.
Args: GENCODE_GTF OUT_CSV"""
import sys
import re
import csv

GENCODE_GTF = sys.argv[1]
OUT_CSV = sys.argv[2]

rows = []
with open(GENCODE_GTF) as f:
    for line in f:
        if line[0] == "#":
            continue
        fields = line.split("\t", 8)
        if fields[2] != "gene":
            continue
        gene_id = re.search(r'gene_id "([^"]+)"', fields[8]).group(1)
        gene_name = re.search(r'gene_name "([^"]+)"', fields[8]).group(1)
        rows.append((gene_id, gene_name))

with open(OUT_CSV, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["gene_id", "gene_name"])
    writer.writerows(rows)

print(f"genes={len(rows)} -> {OUT_CSV}", file=sys.stderr)
