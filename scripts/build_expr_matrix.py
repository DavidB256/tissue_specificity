"""Build a genes x tissues log10(TPM) matrix, restricted to genes with TPM>=1 in
every tissue (a consistent, complete gene set across all 68 tissues)."""
import sys
import pandas as pd
import numpy as np

IN_GCT = sys.argv[1]
OUT_PARQUET = sys.argv[2]

df = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = df.columns[2:]
mat = df[tissues].values
mask_all = (mat >= 1).all(axis=1)
log_mat = np.log10(mat[mask_all])

out = pd.DataFrame(log_mat, columns=tissues)
out.insert(0, "gene", df.loc[mask_all, "Name"].values)
out.to_parquet(OUT_PARQUET, index=False)
print(f"genes={mask_all.sum()} (of {len(df)}) tissues={len(tissues)}", file=sys.stderr)
