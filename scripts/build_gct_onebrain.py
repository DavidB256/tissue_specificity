"""Like build_gct_simplified.py, but collapses ONLY the Brain_* tissue columns into a
single "Brain" column (summing raw TPM, pre floor_tpm) -- every other tissue is left
untouched, unlike build_gct_simplified.py which groups all tissues by first word. Output
is a standard 2-line-header GCT (matches the input format)."""
import sys
import gzip
import pandas as pd

IN_GCT = sys.argv[1]
OUT_GCT = sys.argv[2]

df = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(df.columns[2:])

brain = [t for t in tissues if t.split("_")[0] == "Brain"]
other = [t for t in tissues if t not in brain]

out = df[["Name", "Description"]].copy()
for t in other:
    out[t] = df[t]
out["Brain"] = df[brain].sum(axis=1)

new_tissues = sorted(other) + ["Brain"]
out = out[["Name", "Description"] + new_tissues]

print(f"Brain <- {brain}", file=sys.stderr)
print(f"tissues {len(tissues)} -> {len(new_tissues)}", file=sys.stderr)

with gzip.open(OUT_GCT, "wt") as f:
    f.write("#1.2\n")
    f.write(f"{len(out)}\t{len(new_tissues)}\n")
    out.to_csv(f, sep="\t", index=False)
