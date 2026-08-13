"""Unweighted tau computed on tissue subsets = all non-brain tissues + exactly one
brain tissue, for each of the 13 GTEx Brain_* tissues in turn. Isolates how much
brain's 13-way subregion redundancy inflates unweighted tau, without any GLS/ridge
weighting. Output: one tau column per retained brain tissue, their per-gene mean
(tau_neural_downsampled_mean), and n_nonbrain_tpm1 (count of non-brain tissues with
TPM>=1 per gene, for downstream filtering)."""
import sys
import numpy as np
import pandas as pd
from gls_tau import weighted_tau, floor_tpm

IN_GCT = sys.argv[1]
OUT_CSV = sys.argv[2]

df = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(df.columns[2:])
genes = df["Name"].values
mat = floor_tpm(df[tissues].values.astype(float))

expressed = mat.max(axis=1) > 0
print(f"dropping {(~expressed).sum()}/{len(genes)} genes with no tissue TPM>=1", file=sys.stderr)
genes, mat = genes[expressed], mat[expressed]

brain_idx = [i for i, t in enumerate(tissues) if "Brain" in t]
nonbrain_idx = [i for i in range(len(tissues)) if i not in brain_idx]
brain_tissues = [tissues[i] for i in brain_idx]

out = {
    "gene": genes,
    "n_nonbrain_tpm1": (mat[:, nonbrain_idx] > 0).sum(axis=1),
}
tau_cols = []
for bi, bname in zip(brain_idx, brain_tissues):
    idx = nonbrain_idx + [bi]
    col = f"tau_{bname}"
    out[col] = weighted_tau(mat[:, idx], np.ones(len(idx)))
    tau_cols.append(col)
    print(f"{bname}: done", file=sys.stderr)

out_df = pd.DataFrame(out)
out_df["tau_neural_downsampled_mean"] = out_df[tau_cols].mean(axis=1)
out_df.to_csv(OUT_CSV, index=False)
print(f"genes={len(out_df)} brain_tissues={len(brain_tissues)}", file=sys.stderr)
