"""Bootstrap variance of tau_weighted_gls_ridge using the elementwise-squared tissue
correlation matrix R**2 (see compute_weighted_tau_squared.py) instead of R. Same
resampling schemes as compute_tau_bootstrap_variance.py (tissue / cluster10 / cluster20);
only the ridge tau type is computed here since tau_unweighted never touches R -- reuse
the existing tau_bootstrap_variance_unweighted_{mode}.csv files as-is for that side of the
comparison. Column names match compute_tau_bootstrap_variance.py's ridge output
(tau_ridge_{mode}_var/_mean) so plot_tau_bootstrap_variance.py can be reused unmodified."""
import sys
import numpy as np
import pandas as pd
from gls_tau import tissue_corr_expressed_only, ridge_lambda_search, weighted_tau, floor_tpm

IN_GCT = sys.argv[1]
IN_CLUSTERS = sys.argv[2]
RESAMPLE_MODE = sys.argv[3]  # tissue | cluster10 | cluster20
N_REPS = int(sys.argv[4])
OUT_CSV = sys.argv[5]

SEED = 0

df = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(df.columns[2:])
genes = df["Name"].values
mat = floor_tpm(df[tissues].values.astype(float))
n_tissues, n_genes = len(tissues), mat.shape[0]

rng = np.random.default_rng(SEED)

if RESAMPLE_MODE.startswith("cluster"):
    clusters = pd.read_csv(IN_CLUSTERS)
    tissue_to_idx = {t: i for i, t in enumerate(tissues)}
    blocks = [np.array([tissue_to_idx[t] for t in g["tissue"]])
              for _, g in clusters.groupby("block")]

reps = np.full((N_REPS, n_genes), np.nan)
for r in range(N_REPS):
    if RESAMPLE_MODE == "tissue":
        idx = rng.integers(0, n_tissues, size=n_tissues)
    else:
        block_draw = rng.integers(0, len(blocks), size=len(blocks))
        idx = np.concatenate([blocks[b] for b in block_draw])

    mat_r = mat[:, idx]
    R = tissue_corr_expressed_only(mat_r)
    assert (R >= 0).all(), "squaring flips sign of negative correlations -- not valid here"
    _, w = ridge_lambda_search(R ** 2)
    tau_r = weighted_tau(mat_r, w)
    reps[r] = np.where(np.isfinite(tau_r), tau_r, np.nan)

    if (r + 1) % max(1, N_REPS // 10) == 0:
        print(f"ridge_squared/{RESAMPLE_MODE}: rep {r + 1}/{N_REPS}", file=sys.stderr)

out = pd.DataFrame({
    "gene": genes,
    f"tau_ridge_{RESAMPLE_MODE}_var": np.nanvar(reps, axis=0, ddof=1),
    f"tau_ridge_{RESAMPLE_MODE}_mean": np.nanmean(reps, axis=0),
})
out.to_csv(OUT_CSV, index=False)
print(f"genes={n_genes} reps={N_REPS} tau_type=ridge_squared mode={RESAMPLE_MODE}", file=sys.stderr)
