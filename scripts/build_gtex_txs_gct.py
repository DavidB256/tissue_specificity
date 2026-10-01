"""Transpose the gtex_txs tissue x transcript median-TPM matrix
(build_gtex_txs_median_tpm.py) into a transcript x tissue GCT, parallel to
build_baron_gct.py. Description = gene_id (from build_chess_transcript_meta.py),
unlike Baron's Name==Description -- needed downstream to group transcripts by gene
(compute_multitranscript_variance.py)."""
import sys
import gzip
import pandas as pd

IN_CSV = sys.argv[1]
META_CSV = sys.argv[2]
OUT_GCT = sys.argv[3]

df = pd.read_csv(IN_CSV, index_col=0)  # tissue x transcript
mat = df.T  # transcript x tissue
tissues = list(mat.columns)

gene_id = pd.read_csv(META_CSV).set_index("transcript_id")["gene_id"]

out = pd.DataFrame({"Name": mat.index, "Description": gene_id.reindex(mat.index).values})
for t in tissues:
    out[t] = mat[t].values

with gzip.open(OUT_GCT, "wt") as f:
    f.write("#1.2\n")
    f.write(f"{len(out)}\t{len(tissues)}\n")
    out.to_csv(f, sep="\t", index=False)

print(f"transcripts={len(out)} tissues={len(tissues)} -> {OUT_GCT}", file=sys.stderr)
