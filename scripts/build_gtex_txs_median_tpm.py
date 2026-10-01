"""Combine per-tissue CHESS transcript count tables (GTEx salmon quant, raw/estimated
counts -- NOT TPM despite the "tx_counts.txt" filename, confirmed by column sums ~8.5M-
53M rather than the constant 1e6 a true TPM column would sum to) into one tissue x
transcript median-TPM matrix. Per sample (column), counts -> TPM via CHESS transcript
length (rate=count/length_kb, TPM=rate*1e6/sum(rate) -- same convention as
build_baron_celltype_tpm.py), restricted to transcripts with a known CHESS length
(build_chess_transcript_meta.py; both the TPM numerator and the sum(rate) denominator
drop transcripts missing a length, consistent with Baron's gene-length handling).
Per tissue, the reported value is the median (not mean -- requested explicitly, unlike
Baron's per-cell-type mean) TPM across all patient samples for that tissue."""
import os
import sys
import numpy as np
import pandas as pd

META_CSV = sys.argv[1]
OUT_CSV = sys.argv[2]
COUNTS_FILES = sys.argv[3:]

meta = pd.read_csv(META_CSV).dropna(subset=["length_bp"]).set_index("transcript_id")
length_bp = meta["length_bp"]

tissue_medians = {}
for path in COUNTS_FILES:
    tissue = os.path.basename(os.path.dirname(path))
    counts = pd.read_csv(path, sep="\t", index_col=0)
    n_rows = len(counts.index)
    keep = length_bp.index.intersection(counts.index)
    counts = counts.loc[keep]
    length_kb = (length_bp.loc[keep] / 1000.0).to_numpy()

    rate = counts.to_numpy(dtype=float) / length_kb[:, None]
    total_rate = rate.sum(axis=0)
    tpm = rate * (1e6 / total_rate[None, :])

    median_tpm = np.median(tpm, axis=1)
    tissue_medians[tissue] = pd.Series(median_tpm, index=keep)
    print(f"{tissue}: {counts.shape[1]} samples, {len(keep)}/{n_rows} transcripts with known length",
          file=sys.stderr)

out = pd.DataFrame(tissue_medians).T  # tissue x transcript
out.index.name = "tissue"
out.to_csv(OUT_CSV)
print(f"wrote {out.shape[0]} tissues x {out.shape[1]} transcripts -> {OUT_CSV}", file=sys.stderr)
