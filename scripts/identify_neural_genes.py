"""Genes whose greatest median expression (across all GTEx tissues) is in a Brain_*
tissue -- the case-study set for assessing tau weighting, since the 13 GTEx brain
subregions are highly correlated. Ties for the max go to the first occurrence."""
import sys
import numpy as np
import pandas as pd

IN_GCT = sys.argv[1]
OUT_CSV = sys.argv[2]

df = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = np.array(df.columns[2:])
genes = df["Name"].values
mat = df[list(tissues)].values.astype(float)

max_tpm = mat.max(axis=1)
max_tissue = tissues[mat.argmax(axis=1)]
is_neural = np.char.startswith(max_tissue.astype(str), "Brain_")

out = pd.DataFrame({"gene": genes, "max_tissue": max_tissue, "max_tpm": max_tpm, "is_neural": is_neural})
out.to_csv(OUT_CSV, index=False)
print(f"neural genes (max expression in a Brain_* tissue): {int(is_neural.sum())}/{len(out)}", file=sys.stderr)
