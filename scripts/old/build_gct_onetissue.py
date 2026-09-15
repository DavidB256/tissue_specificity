"""Collapse every tissue column whose first underscore-token equals PREFIX (e.g. "Brain"
matches all Brain_* columns; "Liver" matches Liver, Liver_Hepatocyte, Liver_Mixed_Cell,
Liver_Portal_Tract) into a single column named PREFIX -- every other tissue is left
untouched. Uses raw TPM (not sum), taken before floor_tpm (TPM<1->0, applied downstream
in compute_weighted_tau.py etc), so the collapsed column stays on the same scale as
every other column, which is itself a single un-summed per-tissue median TPM. AGG_METHOD
selects mean (dilutes toward the cross-tissue average) or max (keeps the peak value,
avoiding dilution for genes specific to only a subset of the matched tissues) across the
matched columns. Output is a standard 2-line-header GCT (matches the input format)."""
import sys
import gzip
import pandas as pd

IN_GCT = sys.argv[1]
OUT_GCT = sys.argv[2]
PREFIX = sys.argv[3]
AGG_METHOD = sys.argv[4] if len(sys.argv) > 4 else "mean"
assert AGG_METHOD in ("mean", "max"), f"AGG_METHOD must be mean or max, got {AGG_METHOD!r}"

df = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(df.columns[2:])

matched = [t for t in tissues if t.split("_")[0] == PREFIX]
other = [t for t in tissues if t not in matched]

out = df[["Name", "Description"]].copy()
for t in other:
    out[t] = df[t]
out[PREFIX] = getattr(df[matched], AGG_METHOD)(axis=1)

new_tissues = sorted(other) + [PREFIX]
out = out[["Name", "Description"] + new_tissues]

print(f"{PREFIX} <- {matched}", file=sys.stderr)
print(f"tissues {len(tissues)} -> {len(new_tissues)}", file=sys.stderr)

with gzip.open(OUT_GCT, "wt") as f:
    f.write("#1.2\n")
    f.write(f"{len(out)}\t{len(new_tissues)}\n")
    out.to_csv(f, sep="\t", index=False)
