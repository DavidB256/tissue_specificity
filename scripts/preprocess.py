"""Parse GTEx median-TPM GCT, filter TPM>=1 per tissue, log10, write long table."""
import sys
import pandas as pd
import numpy as np

IN_GCT = sys.argv[1]
OUT_PARQUET = sys.argv[2]

df = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = df.columns[2:]

records = []
for t in tissues:
    v = df[t]
    v = v[v >= 1]
    lv = np.log10(v)
    records.append(pd.DataFrame({"tissue": t, "gene": df.loc[lv.index, "Name"].values, "log10_tpm": lv.values}))

long_df = pd.concat(records, ignore_index=True)
long_df.to_parquet(OUT_PARQUET, index=False)
print(f"tissues={len(tissues)} rows={len(long_df)}", file=sys.stderr)
