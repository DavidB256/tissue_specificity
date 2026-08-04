"""Build a long-format (tissue=condition, gene=transcript, log10_tpm) parquet from
per-condition salmon quant.sf files -- same schema preprocess.py produces, so
fit_one.py works on it unchanged with condition names standing in for tissues."""
import sys
import pandas as pd
import numpy as np

args = sys.argv[1:-1]
OUT_PARQUET = sys.argv[-1]
pairs = list(zip(args[0::2], args[1::2]))

records = []
for condition, quantsf in pairs:
    df = pd.read_csv(quantsf, sep="\t")
    tpm = df["TPM"]
    tpm = tpm[tpm >= 1]
    lv = np.log10(tpm)
    records.append(pd.DataFrame({
        "tissue": condition,
        "gene": df.loc[lv.index, "Name"].values,
        "log10_tpm": lv.values,
    }))

long_df = pd.concat(records, ignore_index=True)
long_df.to_parquet(OUT_PARQUET, index=False)
print(long_df.groupby("tissue").size(), file=sys.stderr)
