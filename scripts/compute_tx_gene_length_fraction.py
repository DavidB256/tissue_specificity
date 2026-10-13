"""Per transcript: tx_length (sum of its own exon lengths) / gene_exonic_length
(union of all exon intervals across all the gene's transcripts, merged).
From a GENCODE GTF."""
import sys
import re
from collections import defaultdict
import pandas as pd

GTF = sys.argv[1]
OUT_CSV = sys.argv[2]

transcript_to_gene = {}
transcript_exons = defaultdict(list)  # transcript_id -> [(start, end)]
gene_exons = defaultdict(list)  # gene_id -> [(start, end)]

with open(GTF) as f:
    for line in f:
        if line[0] == "#":
            continue
        fields = line.rstrip("\n").split("\t", 8)
        if fields[2] != "exon":
            continue
        attrs = fields[8]
        gene_id = re.search(r'gene_id "([^"]+)"', attrs).group(1)
        transcript_id = re.search(r'transcript_id "([^"]+)"', attrs).group(1)
        start, end = int(fields[3]), int(fields[4])

        transcript_to_gene[transcript_id] = gene_id
        transcript_exons[transcript_id].append((start, end))
        gene_exons[gene_id].append((start, end))


def merged_length(intervals):
    intervals = sorted(intervals)
    total = 0
    cur_start, cur_end = intervals[0]
    for start, end in intervals[1:]:
        if start > cur_end + 1:
            total += cur_end - cur_start + 1
            cur_start, cur_end = start, end
        else:
            cur_end = max(cur_end, end)
    total += cur_end - cur_start + 1
    return total


gene_exonic_length = {gene_id: merged_length(intervals) for gene_id, intervals in gene_exons.items()}

rows = []
for transcript_id, exons in transcript_exons.items():
    gene_id = transcript_to_gene[transcript_id]
    tx_length = merged_length(exons)
    gene_length = gene_exonic_length[gene_id]
    rows.append({
        "transcript_id": transcript_id,
        "gene_id": gene_id,
        "tx_length": tx_length,
        "gene_exonic_length": gene_length,
        "length_fraction": tx_length / gene_length,
    })

out = pd.DataFrame(rows)
out.to_csv(OUT_CSV, index=False)

print(f"n_transcripts={len(out)} n_genes={out['gene_id'].nunique()} -> {OUT_CSV}", file=sys.stderr)
