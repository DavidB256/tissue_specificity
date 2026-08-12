"""Rank-stability counterpart to compute_tau_bootstrap_variance.py: instead of per-gene
variance across resamples, measure per-rep Spearman rho between the resampled tau
ranking and the canonical (full-data) tau ranking from weighted_tau.csv. Shows whether
GLS/ridge weighting preserves gene specificity RANKINGS under resampling, vs just
reducing per-gene variance. Same three resampling schemes as compute_tau_bootstrap_variance.py
(tissue / cluster10 / cluster20) and same fixed SEED, so unweighted/ridge jobs for a
given mode see identical resample draws and are directly comparable rep-for-rep."""
import sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from gls_tau import tissue_corr, ridge_lambda_search, weighted_tau, floor_tpm

IN_GCT = sys.argv[1]
IN_TAU = sys.argv[2]
IN_CLUSTERS = sys.argv[3]
TAU_TYPE = sys.argv[4]       # unweighted | ridge
RESAMPLE_MODE = sys.argv[5]  # tissue | cluster10 | cluster20
N_REPS = int(sys.argv[6])
OUT_CSV = sys.argv[7]

SEED = 0
TAU_COL = "tau_unweighted" if TAU_TYPE == "unweighted" else "tau_weighted_gls_ridge"

df = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(df.columns[2:])
genes = df["Name"].values
mat = floor_tpm(df[tissues].values.astype(float))
n_tissues = len(tissues)

tau_ref = pd.read_csv(IN_TAU).set_index("gene").loc[genes, TAU_COL].values

rng = np.random.default_rng(SEED)

if RESAMPLE_MODE.startswith("cluster"):
    clusters = pd.read_csv(IN_CLUSTERS)
    tissue_to_idx = {t: i for i, t in enumerate(tissues)}
    blocks = [np.array([tissue_to_idx[t] for t in g["tissue"]])
              for _, g in clusters.groupby("block")]

rhos = np.full(N_REPS, np.nan)
n_used = np.full(N_REPS, 0)
for r in range(N_REPS):
    if RESAMPLE_MODE == "tissue":
        idx = rng.integers(0, n_tissues, size=n_tissues)
    else:
        block_draw = rng.integers(0, len(blocks), size=len(blocks))
        idx = np.concatenate([blocks[b] for b in block_draw])

    mat_r = mat[:, idx]
    if TAU_TYPE == "unweighted":
        tau_r = weighted_tau(mat_r, np.ones(len(idx)))
    else:
        R = tissue_corr(mat_r)
        _, w = ridge_lambda_search(R)
        tau_r = weighted_tau(mat_r, w)

    ok = np.isfinite(tau_ref) & np.isfinite(tau_r)
    rho, _ = spearmanr(tau_ref[ok], tau_r[ok])
    rhos[r] = rho
    n_used[r] = ok.sum()

    if (r + 1) % max(1, N_REPS // 10) == 0:
        print(f"{TAU_TYPE}/{RESAMPLE_MODE}: rep {r + 1}/{N_REPS}", file=sys.stderr)

out = pd.DataFrame({
    "rep": np.arange(N_REPS),
    f"spearman_rho_{TAU_TYPE}_{RESAMPLE_MODE}": rhos,
    "n_genes": n_used,
})
out.to_csv(OUT_CSV, index=False)
print(f"reps={N_REPS} tau_type={TAU_TYPE} mode={RESAMPLE_MODE} "
      f"median_rho={np.nanmedian(rhos):.4f}", file=sys.stderr)
