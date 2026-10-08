"""Mutually exclusive exons (MXE) in a GENCODE GTF, classic definition: pairs
of exons of the same gene that (i) never co-occur together in any single
transcript (isoform) of that gene, AND (ii) share a flanking context -- some
transcript has exon1 with the same immediate upstream/downstream neighbor
exons (by coords) as exon2 has in some other transcript. Condition (ii)
restricts to the canonical upstream->{A or B}->downstream splicing pattern,
excluding incidental non-co-occurring exon pairs (e.g. unrelated alternative
first/last exons of a multi-promoter/multi-TES gene)."""
import sys
import re
from itertools import combinations
from collections import defaultdict
import pandas as pd

GTF = sys.argv[1]
OUT_CSV = sys.argv[2]
CHROM = sys.argv[3] if len(sys.argv) > 3 else None

genes = {}  # gene_id -> (gene_name, chrom, start, end, strand)
transcript_to_gene = {}
transcript_exons = defaultdict(list)  # transcript_id -> [(exon_number:int, (start, end))]
gene_exons = defaultdict(set)  # gene_id -> {(start, end)}
exon_numbers = defaultdict(set)  # (gene_id, start, end) -> {exon_number}

with open(GTF) as f:
    for line in f:
        if line[0] == "#":
            continue
        fields = line.rstrip("\n").split("\t", 8)
        chrom = fields[0]
        if CHROM and chrom != CHROM:
            continue
        feature = fields[2]
        attrs = fields[8]

        if feature == "gene":
            gene_id = re.search(r'gene_id "([^"]+)"', attrs).group(1)
            gene_name = re.search(r'gene_name "([^"]+)"', attrs).group(1)
            genes[gene_id] = (gene_name, chrom, int(fields[3]), int(fields[4]), fields[6])

        elif feature == "exon":
            gene_id = re.search(r'gene_id "([^"]+)"', attrs).group(1)
            transcript_id = re.search(r'transcript_id "([^"]+)"', attrs).group(1)
            exon_number = int(re.search(r'exon_number (\d+)', attrs).group(1))
            start, end = int(fields[3]), int(fields[4])

            transcript_to_gene[transcript_id] = gene_id
            transcript_exons[transcript_id].append((exon_number, (start, end)))
            gene_exons[gene_id].add((start, end))
            exon_numbers[(gene_id, start, end)].add(str(exon_number))

# flanking signature per exon occurrence: (prev_exon_coords_or_None, next_exon_coords_or_None)
flank_sig = defaultdict(set)  # (gene_id, (start, end)) -> {(prev, next)}
for transcript_id, exon_list in transcript_exons.items():
    gene_id = transcript_to_gene[transcript_id]
    ordered = [e for _, e in sorted(exon_list)]
    for i, exon in enumerate(ordered):
        prev_exon = ordered[i - 1] if i > 0 else None
        next_exon = ordered[i + 1] if i < len(ordered) - 1 else None
        flank_sig[(gene_id, exon)].add((prev_exon, next_exon))

rows = []
for gene_id, all_exons in gene_exons.items():
    if gene_id not in genes or len(all_exons) < 2:
        continue
    gene_name, chrom, gstart, gend, strand = genes[gene_id]

    transcripts = [t for t, g in transcript_to_gene.items() if g == gene_id]
    cooccur = set()
    for t in transcripts:
        exons_t = sorted(e for _, e in transcript_exons[t])
        for e1, e2 in combinations(exons_t, 2):
            cooccur.add((e1, e2))

    for e1, e2 in combinations(sorted(all_exons), 2):
        if (e1, e2) in cooccur:
            continue
        if not (flank_sig[(gene_id, e1)] & flank_sig[(gene_id, e2)]):
            continue
        rows.append({
            "gene_id": gene_id,
            "gene_name": gene_name,
            "gene_chrom": chrom,
            "gene_start": gstart,
            "gene_end": gend,
            "strand": strand,
            "exon1_number": ",".join(sorted(exon_numbers[(gene_id, *e1)], key=int)),
            "exon1_start": e1[0],
            "exon1_end": e1[1],
            "exon2_number": ",".join(sorted(exon_numbers[(gene_id, *e2)], key=int)),
            "exon2_start": e2[0],
            "exon2_end": e2[1],
        })

out = pd.DataFrame(rows)
out.to_csv(OUT_CSV, index=False)

print(f"chrom={CHROM or 'all'} genes_considered={len(gene_exons)} "
      f"MXE_pairs={len(out)} -> {OUT_CSV}", file=sys.stderr)
