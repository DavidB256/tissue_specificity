"""Build a genes x tissues log10(TPM) matrix using ALL genes (no filtering):
TPM>=1 -> log10(TPM); TPM<1 (including 0) -> 0 in log-space."""
import sys
import pandas as pd
import numpy as np

IN_GCT = sys.argv[1]
OUT_PARQUET = sys.argv[2]

df = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = df.columns[2:]
mat = df[tissues].values
log_mat = np.where(mat >= 1, np.log10(np.maximum(mat, 1e-300)), 0.0)

out = pd.DataFrame(log_mat, columns=tissues)
out.insert(0, "gene", df["Name"].values)
out.to_parquet(OUT_PARQUET, index=False)
print(f"genes={len(df)} tissues={len(tissues)} (all genes, TPM<1 floored to 0)", file=sys.stderr)
