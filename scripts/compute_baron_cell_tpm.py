"""Per-cell (not per-cell-type) protein-coding TPM matrix for the 4 Baron et al.
(GSE84133) human pancreas scRNA-seq samples -- same per-cell TPM computation as
build_baron_celltype_tpm.py (gene lengths from build_baron_gene_lengths.py),
but stops short of that script's final groupby-mean-by-cell-type step, and
additionally restricts to protein-coding genes (GENCODE gene_name match, same
convention as filter_protein_coding.py). Keeps assigned_cluster (Baron's own
cell-type label) and sample columns for reference -- not used to build the
matrix itself."""
import sys
import re
import numpy as np
import pandas as pd

LENGTHS_CSV = sys.argv[-3]
GENCODE_GTF = sys.argv[-2]
OUT_CSV = sys.argv[-1]
SAMPLE_CSVS = sys.argv[1:-3]

lengths = pd.read_csv(LENGTHS_CSV)
lengths = lengths[lengths["source"] != "neither"].set_index("gene")["length_bp"]
known_genes = lengths.index

pc_genes = set()
with open(GENCODE_GTF) as f:
    for line in f:
        if line[0] == "#":
            continue
        fields = line.split("\t", 8)
        if fields[2] != "gene" or 'gene_type "protein_coding"' not in fields[8]:
            continue
        pc_genes.add(re.search(r'gene_name "([^"]+)"', fields[8]).group(1))

frames = []
sample_labels = []
for path in SAMPLE_CSVS:
    df = pd.read_csv(path, index_col=0)
    keep_genes = [g for g in known_genes if g in df.columns and g in pc_genes]
    frames.append(df[["assigned_cluster"] + keep_genes])
    sample_labels.extend([path] * len(df))
    print(f"{path}: {len(df)} cells, {len(keep_genes)} protein-coding genes kept", file=sys.stderr)

combined = pd.concat(frames, axis=0, ignore_index=True)
combined["sample"] = sample_labels
gene_cols = [c for c in combined.columns if c not in ("assigned_cluster", "sample")]
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
tpm_df.insert(0, "sample", combined["sample"].values)
tpm_df.insert(0, "assigned_cluster", combined["assigned_cluster"].values)

tpm_df.to_csv(OUT_CSV, index=False)
print(f"wrote {len(tpm_df)} cells x {len(gene_cols)} genes -> {OUT_CSV}", file=sys.stderr)
