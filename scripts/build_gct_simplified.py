"""Collapse fine-grained GTEx sub-tissue columns into their parent tissue by summing
raw TPM (pre floor_tpm) across all columns whose name is <parent_tissue>_<suffix>, where
<parent_tissue> is itself an existing tissue column (e.g. Liver_Hepatocyte/Mixed_Cell/
Portal_Tract -> Liver; Stomach_Mucosa/Muscularis/Mixed_Cell -> Stomach; Pancreas_Acini/
Islets/Mixed_Cell -> Pancreas; Colon_Transverse_Mixed_Cell/Mucosa/Muscularis ->
Colon_Transverse; Small_Intestine_Terminal_Ileum_Lymphoid_Aggregate/Mixed_Cell ->
Small_Intestine_Terminal_Ileum). Summing raw TPM prior to the TPM<1 floor lets low-TPM
subtype signal "rescue" a gene above the floor once merged into its parent. Tissues with
no matching parent column (e.g. Brain_*, Artery_*, Skin_*, Adipose_*) are left unchanged
-- their names don't syntactically decompose into <existing_tissue>_<suffix>. Output is a
standard 2-line-header GCT (matches the input format) so it's a drop-in input for the
rest of the pipeline."""
import sys
import gzip
import pandas as pd

IN_GCT = sys.argv[1]
OUT_GCT = sys.argv[2]

df = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(df.columns[2:])
colset = set(tissues)


def find_base(col):
    parts = col.split("_")
    for i in range(len(parts) - 1, 0, -1):
        prefix = "_".join(parts[:i])
        if prefix in colset:
            return prefix
    return col


groups = {}
for t in tissues:
    groups.setdefault(find_base(t), []).append(t)

new_tissues = sorted(groups)
out = df[["Name", "Description"]].copy()
for base in new_tissues:
    out[base] = df[groups[base]].sum(axis=1)

for base, members in groups.items():
    if len(members) > 1:
        print(f"{base} <- {members}", file=sys.stderr)
print(f"tissues {len(tissues)} -> {len(new_tissues)}", file=sys.stderr)

with gzip.open(OUT_GCT, "wt") as f:
    f.write("#1.2\n")
    f.write(f"{len(out)}\t{len(new_tissues)}\n")
    out.to_csv(f, sep="\t", index=False)
