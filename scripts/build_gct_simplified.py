"""Collapse GTEx tissue columns by their first underscore-delimited word (e.g. all
Brain_* -> Brain, all Artery_* -> Artery, Liver/Liver_Hepatocyte/Liver_Mixed_Cell/
Liver_Portal_Tract -> Liver), summing raw TPM (pre floor_tpm) across each group. Summing
prior to the TPM<1 floor lets low-TPM subtype signal "rescue" a gene above the floor once
merged. Groups of size 1 keep their original full tissue name (no truncation when nothing
is being merged). Output is a standard 2-line-header GCT (matches the input format) so
it's a drop-in input for the rest of the pipeline."""
import sys
import gzip
import pandas as pd

IN_GCT = sys.argv[1]
OUT_GCT = sys.argv[2]

df = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(df.columns[2:])

groups = {}
for t in tissues:
    groups.setdefault(t.split("_")[0], []).append(t)

out = df[["Name", "Description"]].copy()
new_tissues = []
for key in sorted(groups):
    members = groups[key]
    name = key if len(members) > 1 else members[0]
    new_tissues.append(name)
    out[name] = df[members].sum(axis=1)

for base, members in groups.items():
    if len(members) > 1:
        print(f"{base} <- {members}", file=sys.stderr)
print(f"tissues {len(tissues)} -> {len(new_tissues)}", file=sys.stderr)

with gzip.open(OUT_GCT, "wt") as f:
    f.write("#1.2\n")
    f.write(f"{len(out)}\t{len(new_tissues)}\n")
    out.to_csv(f, sep="\t", index=False)
