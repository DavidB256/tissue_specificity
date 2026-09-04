"""Build a gene -> length(bp) table for Baron et al. (GSE84133) pancreas scRNA-seq
genes, for TPM conversion. Length = summed exon length of the MANE_Select
transcript if the gene is in MANE; else the longest transcript (by summed exon
bp, not genomic span) in GENCODE; else "neither" (length left blank, gene
dropped downstream)."""
import re
import sys

import pandas as pd

MANE_GTF = sys.argv[1]
GENCODE_GTF = sys.argv[2]
BARON_CSV = sys.argv[3]
OUT_CSV = sys.argv[4]

_ATTR_RE = {}


def attr(attrs, key):
    pat = _ATTR_RE.get(key)
    if pat is None:
        pat = re.compile(key + r' "([^"]+)"')
        _ATTR_RE[key] = pat
    m = pat.search(attrs)
    return m.group(1) if m else None


def exon_lengths_by_transcript(gtf_path, require_mane_select=False):
    """Return {transcript_id: summed_exon_len}, {transcript_id: gene_name}."""
    tx_len, tx_gene = {}, {}
    n = 0
    with open(gtf_path) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[2] != "exon":
                continue
            attrs = fields[8]
            if require_mane_select and 'tag "MANE_Select"' not in attrs:
                continue
            tid = attr(attrs, "transcript_id")
            gname = attr(attrs, "gene_name")
            length = int(fields[4]) - int(fields[3]) + 1
            tx_len[tid] = tx_len.get(tid, 0) + length
            tx_gene[tid] = gname
            n += 1
            if n % 500_000 == 0:
                print(f"  ...{n} exon lines ({gtf_path})", file=sys.stderr)
    return tx_len, tx_gene


def longest_per_gene(tx_len, tx_gene):
    gene_len = {}
    for tid, length in tx_len.items():
        gname = tx_gene[tid]
        if gname not in gene_len or length > gene_len[gname]:
            gene_len[gname] = length
    return gene_len


print("parsing MANE (MANE_Select transcripts only)...", file=sys.stderr)
mane_tx_len, mane_tx_gene = exon_lengths_by_transcript(MANE_GTF, require_mane_select=True)
mane_gene_len = longest_per_gene(mane_tx_len, mane_tx_gene)
print(f"MANE_Select genes: {len(mane_gene_len)}", file=sys.stderr)

print("parsing GENCODE (all transcripts, longest per gene)...", file=sys.stderr)
gc_tx_len, gc_tx_gene = exon_lengths_by_transcript(GENCODE_GTF, require_mane_select=False)
gencode_gene_len = longest_per_gene(gc_tx_len, gc_tx_gene)
print(f"GENCODE genes: {len(gencode_gene_len)}", file=sys.stderr)

baron_genes = list(pd.read_csv(BARON_CSV, index_col=0, nrows=0).columns[2:])

rows = []
n_mane = n_gencode = n_neither = 0
for gene in baron_genes:
    if gene in mane_gene_len:
        rows.append((gene, "MANE", mane_gene_len[gene]))
        n_mane += 1
    elif gene in gencode_gene_len:
        rows.append((gene, "GENCODE", gencode_gene_len[gene]))
        n_gencode += 1
    else:
        rows.append((gene, "neither", None))
        n_neither += 1

out = pd.DataFrame(rows, columns=["gene", "source", "length_bp"])
out.to_csv(OUT_CSV, index=False)

pct_neither = 100 * n_neither / len(baron_genes)
print(
    f"baron genes={len(baron_genes)} MANE={n_mane} GENCODE={n_gencode} "
    f"neither={n_neither} ({pct_neither:.2f}%)",
    file=sys.stderr,
)
if pct_neither > 5:
    print(
        f"WARNING: {pct_neither:.2f}% of Baron genes found in neither MANE nor "
        "GENCODE (>5%) -- these are dropped from the TPM matrix.",
        file=sys.stderr,
    )
