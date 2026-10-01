"""Transcript-level metadata table from the CHESS 3.1.3 annotation GTF: transcript_id,
gene_id, gene_type (both carried on the transcript feature line directly), length_bp
(summed exon length, from the transcript's own exon feature lines -- not a MANE/GENCODE
fallback like build_baron_gene_lengths.py, since CHESS already has one length per
transcript). Used both as the TPM length-normalization denominator
(build_gtex_txs_median_tpm.py) and the protein_coding filter
(filter_protein_coding_transcripts.py)."""
import re
import sys
import pandas as pd

GTF = sys.argv[1]
OUT_CSV = sys.argv[2]

_ATTR_RE = {}


def attr(attrs, key):
    pat = _ATTR_RE.get(key)
    if pat is None:
        pat = re.compile(key + r' "([^"]+)"')
        _ATTR_RE[key] = pat
    m = pat.search(attrs)
    return m.group(1) if m else None


tx_gene, tx_type, tx_len = {}, {}, {}
n_exon = n_tx = 0
with open(GTF) as f:
    for line in f:
        if line[0] == "#":
            continue
        fields = line.rstrip("\n").split("\t", 8)
        attrs = fields[8]
        if fields[2] == "transcript":
            tid = attr(attrs, "transcript_id")
            tx_gene[tid] = attr(attrs, "gene_id")
            tx_type[tid] = attr(attrs, "gene_type")
            n_tx += 1
        elif fields[2] == "exon":
            tid = attr(attrs, "transcript_id")
            length = int(fields[4]) - int(fields[3]) + 1
            tx_len[tid] = tx_len.get(tid, 0) + length
            n_exon += 1
        if (n_exon + n_tx) % 500_000 == 0:
            print(f"  ...{n_tx} transcript / {n_exon} exon lines", file=sys.stderr)

out = pd.DataFrame({
    "transcript_id": list(tx_gene.keys()),
    "gene_id": [tx_gene[t] for t in tx_gene],
    "gene_type": [tx_type[t] for t in tx_gene],
    "length_bp": [tx_len.get(t) for t in tx_gene],
})
out.to_csv(OUT_CSV, index=False)

n_missing_len = out["length_bp"].isna().sum()
print(f"transcripts={len(out)} missing_length={n_missing_len} "
      f"protein_coding={(out['gene_type'] == 'protein_coding').sum()} -> {OUT_CSV}", file=sys.stderr)
