"""Like build_gct_onebrain.py, but collapses ONLY the Liver* tissue columns (Liver,
Liver_Hepatocyte, Liver_Mixed_Cell, Liver_Portal_Tract) into a single "Liver" column
(summing raw TPM, pre floor_tpm) -- every other tissue is left untouched. Output is a
standard 2-line-header GCT (matches the input format)."""
import sys
import gzip
import pandas as pd

IN_GCT = sys.argv[1]
OUT_GCT = sys.argv[2]

df = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(df.columns[2:])

liver = [t for t in tissues if t.split("_")[0] == "Liver"]
other = [t for t in tissues if t not in liver]

out = df[["Name", "Description"]].copy()
for t in other:
    out[t] = df[t]
out["Liver"] = df[liver].sum(axis=1)

new_tissues = sorted(other) + ["Liver"]
out = out[["Name", "Description"] + new_tissues]

print(f"Liver <- {liver}", file=sys.stderr)
print(f"tissues {len(tissues)} -> {len(new_tissues)}", file=sys.stderr)

with gzip.open(OUT_GCT, "wt") as f:
    f.write("#1.2\n")
    f.write(f"{len(out)}\t{len(new_tissues)}\n")
    out.to_csv(f, sep="\t", index=False)
