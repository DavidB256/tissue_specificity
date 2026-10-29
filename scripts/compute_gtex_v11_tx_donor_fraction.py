"""Stream the GTEx v11 RSEM transcript TPM table (see build_gtex_v11_tx_gct.py) and, for
the transcripts in the given GCT (e.g. the protein-coding subset), write the fraction of
samples (donors) per tissue with TPM > 0 -> transcript x tissue CSV."""
import sys
import numpy as np
import pandas as pd

CHUNKSIZE = 4000

IN_TPM_TXT_GZ = sys.argv[1]
IN_SAMPLE_TISSUE_CSV = sys.argv[2]
IN_GCT = sys.argv[3]
OUT_CSV = sys.argv[4]

keep = set(pd.read_csv(IN_GCT, sep="\t", skiprows=2, usecols=["Name"])["Name"])
sample_tissue = pd.read_csv(IN_SAMPLE_TISSUE_CSV).set_index("sample_id")["tissue"]
header = pd.read_csv(IN_TPM_TXT_GZ, sep="\t", nrows=0).columns
sample_cols = [c for c in header[2:] if c in sample_tissue.index]
tissues = sorted(sample_tissue.loc[sample_cols].unique())
tissue_cols = {t: [c for c in sample_cols if sample_tissue[c] == t] for t in tissues}

dtype = {"transcript_id": str}
dtype.update({c: np.float32 for c in sample_cols})
out = []
n = 0
for chunk in pd.read_csv(IN_TPM_TXT_GZ, sep="\t", usecols=["transcript_id"] + sample_cols,
                         dtype=dtype, chunksize=CHUNKSIZE):
    chunk = chunk[chunk["transcript_id"].isin(keep)]
    res = pd.DataFrame({"transcript_id": chunk["transcript_id"].values})
    for t in tissues:
        res[t] = (chunk[tissue_cols[t]].to_numpy() > 0).mean(axis=1)
    out.append(res)
    n += CHUNKSIZE
    print(f"  ...{n} rows scanned", file=sys.stderr)

out = pd.concat(out, ignore_index=True)
out.to_csv(OUT_CSV, index=False)
print(f"transcripts={len(out)} tissues={len(tissues)} -> {OUT_CSV}", file=sys.stderr)
