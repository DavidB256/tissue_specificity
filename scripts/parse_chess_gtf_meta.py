"""transcript_id,gene_id,gene_name,gene_type,chrom,strand for every transcript in
a CHESS gtf (attrs on the "transcript" feature line itself -- no separate "gene"
feature). All ids are "CHS.<gene>.<tx>", no version suffix."""
import sys
import re
import pandas as pd

IN_GTF = sys.argv[1]
OUT_CSV = sys.argv[2]

rows = []
with open(IN_GTF) as f:
    for line in f:
        if line[0] == "#":
            continue
        fields = line.rstrip("\n").split("\t", 8)
        if fields[2] != "transcript":
            continue
        attrs = fields[8]
        transcript_id = re.search(r'transcript_id "([^"]+)"', attrs).group(1)
        gene_id = re.search(r'gene_id "([^"]+)"', attrs).group(1)
        gene_name = re.search(r'gene_name "([^"]+)"', attrs)
        gene_type = re.search(r'gene_type "([^"]+)"', attrs)
        rows.append({
            "transcript_id": transcript_id, "gene_id": gene_id,
            "gene_name": gene_name.group(1) if gene_name else None,
            "gene_type": gene_type.group(1) if gene_type else None,
            "chrom": fields[0], "strand": fields[6],
        })

out = pd.DataFrame(rows)
out.to_csv(OUT_CSV, index=False)

print(f"n_transcripts={len(out)} n_genes={out['gene_id'].nunique()} -> {OUT_CSV}", file=sys.stderr)
