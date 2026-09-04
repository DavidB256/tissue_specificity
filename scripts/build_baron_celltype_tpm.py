"""Combine the 4 Baron et al. (GSE84133) human pancreas scRNA-seq samples into
one cell-type x gene TPM matrix. Per cell, counts are converted to TPM using
gene lengths from build_baron_gene_lengths.py (genes with no MANE/GENCODE
length are dropped from both numerator and the length-normalization
denominator). Per cell type (assigned_cluster), the reported value is the mean
per-cell TPM across all cells of that type, pooled across the 4 samples."""
import sys

import numpy as np
import pandas as pd

SAMPLE_CSVS = sys.argv[1:-2]
LENGTHS_CSV = sys.argv[-2]
OUT_CSV = sys.argv[-1]

lengths = pd.read_csv(LENGTHS_CSV)
lengths = lengths[lengths["source"] != "neither"].set_index("gene")["length_bp"]
known_genes = lengths.index

frames = []
for path in SAMPLE_CSVS:
    df = pd.read_csv(path, index_col=0)
    keep_genes = [g for g in known_genes if g in df.columns]
    frames.append(df[["assigned_cluster"] + keep_genes])
    print(f"{path}: {len(df)} cells, {len(keep_genes)} genes kept", file=sys.stderr)

combined = pd.concat(frames, axis=0, ignore_index=True)
gene_cols = [c for c in combined.columns if c != "assigned_cluster"]
print(f"combined: {len(combined)} cells, {len(gene_cols)} genes", file=sys.stderr)

length_kb = (lengths[gene_cols] / 1000).to_numpy()
counts = combined[gene_cols].to_numpy(dtype=float)
rate = counts / length_kb
total_rate = rate.sum(axis=1)

zero_mask = total_rate == 0
if zero_mask.any():
    print(f"dropping {zero_mask.sum()} cells with zero counts in kept genes", file=sys.stderr)
combined = combined.loc[~zero_mask].reset_index(drop=True)
rate = rate[~zero_mask]
total_rate = total_rate[~zero_mask]

tpm = rate * (1e6 / total_rate[:, None])
tpm_df = pd.DataFrame(tpm, columns=gene_cols)
tpm_df["assigned_cluster"] = combined["assigned_cluster"].values

celltype_tpm = tpm_df.groupby("assigned_cluster").mean()
celltype_tpm.index.name = "cell_type"

counts_per_type = tpm_df["assigned_cluster"].value_counts()
print("cells per cell type:", file=sys.stderr)
print(counts_per_type.to_string(), file=sys.stderr)

celltype_tpm.to_csv(OUT_CSV)
print(f"wrote {celltype_tpm.shape[0]} cell types x {celltype_tpm.shape[1]} genes -> {OUT_CSV}", file=sys.stderr)
