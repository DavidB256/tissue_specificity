"""Subset a GCT to the tissues listed (one per line) in a text file; keeps all rows."""
import sys
import gzip
import pandas as pd

IN_GCT, IN_TISSUES, OUT_GCT = sys.argv[1:4]
keep = [l.strip() for l in open(IN_TISSUES) if l.strip()]
df = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
missing = [t for t in keep if t not in df.columns]
assert not missing, f"tissues not in GCT: {missing}"
df = df[list(df.columns[:2]) + [c for c in df.columns[2:] if c in set(keep)]]
with gzip.open(OUT_GCT, "wt") as f:
    f.write(f"#1.2\n{len(df)}\t{df.shape[1] - 2}\n")
    df.to_csv(f, sep="\t", index=False)
print(f"rows={len(df)} tissues={df.shape[1] - 2} -> {OUT_GCT}", file=sys.stderr)
