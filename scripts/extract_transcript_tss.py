"""Transcript ID -> TSS (transcription start site) from a GENCODE GTF.
TSS = transcript start coord on + strand, transcript end coord on - strand
(GTF coords are already 1-based, strand-agnostic start<end)."""
import sys
import re
import pandas as pd

GTF = sys.argv[1]
OUT_CSV = sys.argv[2]

rows = []
with open(GTF) as f:
    for line in f:
        if line[0] == "#":
            continue
        fields = line.rstrip("\n").split("\t", 8)
        if fields[2] != "transcript":
            continue
        attrs = fields[8]
        transcript_id = re.search(r'transcript_id "([^"]+)"', attrs).group(1)
        chrom, strand = fields[0], fields[6]
        start, end = int(fields[3]), int(fields[4])
        tss = start if strand == "+" else end
        rows.append({"transcript_id": transcript_id, "chrom": chrom, "strand": strand, "tss": tss})

out = pd.DataFrame(rows)
out.to_csv(OUT_CSV, index=False)

print(f"n_transcripts={len(out)} -> {OUT_CSV}", file=sys.stderr)
