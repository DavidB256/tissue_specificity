"""Combine per-tissue, per-sample TPM gct.gz files (many samples per tissue) into one
gene x tissue-median GCT.gz, matching the human GTEx median-TPM format so it can be
fed straight into preprocess.py unchanged."""
import sys
import io
import gzip
import pandas as pd

TISSUES = sys.argv[1].split(",")
IN_FILES = sys.argv[2:-1]
OUT_GCT = sys.argv[-1]
assert len(TISSUES) == len(IN_FILES), f"{len(TISSUES)} tissues vs {len(IN_FILES)} files"

combined = None
for tissue, path in zip(TISSUES, IN_FILES):
    df = pd.read_csv(path, sep="\t", skiprows=2)
    sample_cols = df.columns[2:]
    tissue_df = df[["Name", "Description"]].copy()
    tissue_df[tissue] = df[sample_cols].median(axis=1)
    combined = tissue_df if combined is None else combined.merge(tissue_df, on=["Name", "Description"], how="outer")

buf = io.StringIO()
buf.write("#1.2\n")
buf.write(f"{len(combined)}\t{len(TISSUES)}\n")
combined.to_csv(buf, sep="\t", index=False)
with gzip.open(OUT_GCT, "wt") as f:
    f.write(buf.getvalue())

print(f"genes={len(combined)} tissues={len(TISSUES)}", file=sys.stderr)
