"""Like build_gct_onebrain.py, but collapses ONLY the Kidney_* tissue columns into a
single "Kidney" column (summing raw TPM, pre floor_tpm) -- every other tissue is left
untouched. Output is a standard 2-line-header GCT (matches the input format)."""
import sys
import gzip
import pandas as pd

IN_GCT = sys.argv[1]
OUT_GCT = sys.argv[2]

df = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(df.columns[2:])

kidney = [t for t in tissues if t.split("_")[0] == "Kidney"]
other = [t for t in tissues if t not in kidney]

out = df[["Name", "Description"]].copy()
for t in other:
    out[t] = df[t]
out["Kidney"] = df[kidney].sum(axis=1)

new_tissues = sorted(other) + ["Kidney"]
out = out[["Name", "Description"] + new_tissues]

print(f"Kidney <- {kidney}", file=sys.stderr)
print(f"tissues {len(tissues)} -> {len(new_tissues)}", file=sys.stderr)

with gzip.open(OUT_GCT, "wt") as f:
    f.write("#1.2\n")
    f.write(f"{len(out)}\t{len(new_tissues)}\n")
    out.to_csv(f, sep="\t", index=False)
