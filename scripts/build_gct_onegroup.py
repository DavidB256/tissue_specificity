"""Collapse an explicit, named set of columns (MEMBERS) into a single column named
GROUP_NAME -- every other column is left untouched. Generalizes build_gct_onetissue.py's
underscore-prefix matching (e.g. GTEx's Brain_* subregions) to columns that don't share a
common prefix, e.g. Baron pancreas islet cell types (alpha, beta, gamma, delta -> islets).
Uses raw TPM (not sum), taken before floor_tpm (TPM<1->0, applied downstream in
compute_weighted_tau.py etc). AGG_METHOD selects mean or max across MEMBERS. Output is a
standard 2-line-header GCT (matches the input format)."""
import sys
import gzip
import pandas as pd

IN_GCT = sys.argv[1]
OUT_GCT = sys.argv[2]
GROUP_NAME = sys.argv[3]
AGG_METHOD = sys.argv[4]
MEMBERS = sys.argv[5:]
assert AGG_METHOD in ("mean", "max"), f"AGG_METHOD must be mean or max, got {AGG_METHOD!r}"

df = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(df.columns[2:])

missing = [m for m in MEMBERS if m not in tissues]
assert not missing, f"MEMBERS not found in {IN_GCT}: {missing}"

other = [t for t in tissues if t not in MEMBERS]

out = df[["Name", "Description"]].copy()
for t in other:
    out[t] = df[t]
out[GROUP_NAME] = getattr(df[MEMBERS], AGG_METHOD)(axis=1)

new_tissues = sorted(other) + [GROUP_NAME]
out = out[["Name", "Description"] + new_tissues]

print(f"{GROUP_NAME} <- {MEMBERS}", file=sys.stderr)
print(f"tissues {len(tissues)} -> {len(new_tissues)}", file=sys.stderr)

with gzip.open(OUT_GCT, "wt") as f:
    f.write("#1.2\n")
    f.write(f"{len(out)}\t{len(new_tissues)}\n")
    out.to_csv(f, sep="\t", index=False)
