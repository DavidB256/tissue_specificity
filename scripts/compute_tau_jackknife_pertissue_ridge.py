"""Per-gene, per-tissue jackknife tau_weighted_gls_ridge: for each gene (TPM>=1 in >=2
tissues) and each tissue, recompute tau_weighted_gls_ridge with that tissue dropped, USING
THE FIXED weight_gls_ridge weights from weighted_tau_weights.csv (solved once from the
full 68-tissue R) -- weights are not resolved per jackknife rep, only the dropped
tissue's weight entry is removed. One row per gene, one column per tissue (jackknife tau
with that tissue left out), plus tau_var and tau_range across all tissue-dropped samples
for that gene. See compute_tau_jackknife_pertissue.py (unweighted counterpart)."""
import sys
import numpy as np
import pandas as pd
from gls_tau import weighted_tau, floor_tpm

IN_GCT = sys.argv[1]
IN_WEIGHTS = sys.argv[2]
OUT_CSV = sys.argv[3]

df = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(df.columns[2:])
genes = df["Name"].values
mat = floor_tpm(df[tissues].values.astype(float))
n_tissues = len(tissues)

w_full = pd.read_csv(IN_WEIGHTS).set_index("tissue").loc[tissues, "weight_gls_ridge"].values

keep = (mat > 0).sum(axis=1) >= 2
print(f"keeping {keep.sum()}/{len(genes)} genes with TPM>=1 in >=2 tissues", file=sys.stderr)
genes, mat = genes[keep], mat[keep]

jack = np.full((mat.shape[0], n_tissues), np.nan)
for i, t in enumerate(tissues):
    idx = np.arange(n_tissues) != i
    mat_r = mat[:, idx]
    w = w_full[idx]
    jack[:, i] = weighted_tau(mat_r, w)
    print(f"drop {i + 1}/{n_tissues}: {t}", file=sys.stderr)

out = pd.DataFrame(jack, columns=tissues)
out.insert(0, "gene", genes)
out["tau_var"] = np.nanvar(jack, axis=1)
out["tau_range"] = np.nanmax(jack, axis=1) - np.nanmin(jack, axis=1)
out.to_csv(OUT_CSV, index=False)
print(f"genes={len(out)} tissues={n_tissues}", file=sys.stderr)
