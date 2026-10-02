"""Stream the (33GB uncompressed) GTEx v11 RSEM transcript TPM table in row chunks
(too large to hold in memory at once: ~300k transcripts x ~19.8k samples) and reduce
each chunk straight to a transcript x tissue GCT -- median TPM across each tissue's
samples (build_gtex_v11_sample_tissue_map.py), TPM<1 floored to 0 immediately after
that per-tissue combination (per-tissue median first, then floor -- requested
explicitly for this source, unlike gtex_txs). Each chunk is reduced from ~19.8k sample
columns down to ~68 tissue columns before being kept in memory, so accumulating all
chunks is cheap even though the raw input isn't."""
import sys
import numpy as np
import pandas as pd

CHUNKSIZE = 4000

IN_TPM_TXT_GZ = sys.argv[1]
IN_SAMPLE_TISSUE_CSV = sys.argv[2]
OUT_GCT = sys.argv[3]

sample_tissue = pd.read_csv(IN_SAMPLE_TISSUE_CSV).set_index("sample_id")["tissue"]

header = pd.read_csv(IN_TPM_TXT_GZ, sep="\t", nrows=0).columns
sample_cols = [c for c in header[2:] if c in sample_tissue.index]
unmapped = [c for c in header[2:] if c not in sample_tissue.index]
if unmapped:
    print(f"WARNING: {len(unmapped)} sample columns have no tissue mapping, dropped: "
          f"{unmapped[:5]}{'...' if len(unmapped) > 5 else ''}", file=sys.stderr)

tissues = sorted(sample_tissue.loc[sample_cols].unique())
tissue_cols = {t: [c for c in sample_cols if sample_tissue[c] == t] for t in tissues}
print(f"sample columns: {len(sample_cols)} mapped / {len(header) - 2} total, "
      f"tissues={len(tissues)}", file=sys.stderr)

dtype = {"transcript_id": str, "gene_id": str}
dtype.update({c: np.float32 for c in sample_cols})

chunks = []
n_rows = 0
reader = pd.read_csv(IN_TPM_TXT_GZ, sep="\t", usecols=["transcript_id", "gene_id"] + sample_cols,
                      dtype=dtype, chunksize=CHUNKSIZE)
for chunk in reader:
    out_chunk = pd.DataFrame({"Name": chunk["transcript_id"].values, "Description": chunk["gene_id"].values})
    for t in tissues:
        med = chunk[tissue_cols[t]].median(axis=1).to_numpy()
        med[med < 1] = 0.0
        out_chunk[t] = med
    chunks.append(out_chunk)
    n_rows += len(chunk)
    print(f"  ...{n_rows} transcripts processed", file=sys.stderr)

out = pd.concat(chunks, ignore_index=True)

import gzip
with gzip.open(OUT_GCT, "wt") as f:
    f.write("#1.2\n")
    f.write(f"{len(out)}\t{len(tissues)}\n")
    out.to_csv(f, sep="\t", index=False)

print(f"transcripts={len(out)} tissues={len(tissues)} -> {OUT_GCT}", file=sys.stderr)
