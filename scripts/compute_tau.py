"""Compute tau (tissue-specificity index, Yanai et al. 2005) per gene from a GTEx
median-TPM GCT. tau = sum_i(1 - x_i/max_i(x_i)) / (n-1), on linear (non-log) TPM --
tau assumes ratios reflect expression scale directly, so it is computed on raw TPM,
not log-transformed values."""
import sys
import pandas as pd
import numpy as np

IN_GCT = sys.argv[1]
OUT_CSV = sys.argv[2]

df = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = df.columns[2:]
mat = df[tissues].values.astype(float)
n = mat.shape[1]

max_i = mat.max(axis=1)
mean_tpm = mat.mean(axis=1)
with np.errstate(invalid="ignore", divide="ignore"):
    xhat = mat / max_i[:, None]
    tau = np.where(max_i > 0, (1 - xhat).sum(axis=1) / (n - 1), np.nan)

out = pd.DataFrame({"gene": df["Name"].values, "tau": tau, "mean_tpm": mean_tpm, "max_tpm": max_i})
out.to_csv(OUT_CSV, index=False)
print(f"genes={len(out)} tissues={n} nan_tau(all-zero genes)={out['tau'].isna().sum()}", file=sys.stderr)
