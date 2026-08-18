"""Per-gene, per-tissue jackknife tau_unweighted: for each gene (TPM>=1 in >=2 tissues)
and each tissue, recompute tau_unweighted with that tissue dropped. One row per gene, one
column per tissue (jackknife tau with that tissue left out), plus tau_var and tau_range
across all tissue-dropped samples for that gene."""
import sys
import numpy as np
import pandas as pd
from gls_tau import weighted_tau, floor_tpm

IN_GCT = sys.argv[1]
OUT_CSV = sys.argv[2]

df = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(df.columns[2:])
genes = df["Name"].values
mat = floor_tpm(df[tissues].values.astype(float))
n_tissues = len(tissues)

keep = (mat > 0).sum(axis=1) >= 2
print(f"keeping {keep.sum()}/{len(genes)} genes with TPM>=1 in >=2 tissues", file=sys.stderr)
genes, mat = genes[keep], mat[keep]

jack = np.full((mat.shape[0], n_tissues), np.nan)
for i, t in enumerate(tissues):
    idx = np.arange(n_tissues) != i
    jack[:, i] = weighted_tau(mat[:, idx], np.ones(n_tissues - 1))
    print(f"drop {i + 1}/{n_tissues}: {t}", file=sys.stderr)

out = pd.DataFrame(jack, columns=tissues)
out.insert(0, "gene", genes)
out["tau_var"] = np.nanvar(jack, axis=1)
out["tau_range"] = np.nanmax(jack, axis=1) - np.nanmin(jack, axis=1)
out.to_csv(OUT_CSV, index=False)
print(f"genes={len(out)} tissues={n_tissues}", file=sys.stderr)
