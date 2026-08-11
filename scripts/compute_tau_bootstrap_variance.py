"""Bootstrap variance of tau_unweighted / tau_weighted_gls_ridge per gene, under two
resampling schemes:
  tissue  -- resample the 68 GTEx tissues with replacement (standard tissue bootstrap).
  cluster -- resample 10 top-level tissue blocks with replacement (cut at t=10, i.e.
             the first 9 splits of the precomputed hierarchical clustering dendrogram
             -- see compute_tissue_correlation.py); each drawn block contributes all
             its member tissues. Tests whether GLS/ridge weighting is specifically
             robust to whole correlated blocks (e.g. all 13 Brain_* tissues) being
             duplicated/dropped together, vs single tissues under plain tissue
             resampling.

Seed is fixed (not a CLI arg) so tissue/cluster draws are identical across the
unweighted and ridge jobs for the same mode -- keeps the two tau types directly
comparable per rep even though they run as separate parallel snakemake jobs."""
import sys
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster
from gls_tau import tissue_corr, ridge_lambda_search, weighted_tau, floor_tpm

IN_GCT = sys.argv[1]
IN_LINKAGE = sys.argv[2]
TAU_TYPE = sys.argv[3]       # unweighted | ridge
RESAMPLE_MODE = sys.argv[4]  # tissue | cluster
N_REPS = int(sys.argv[5])
OUT_CSV = sys.argv[6]

SEED = 0
N_BLOCKS = 10

df = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(df.columns[2:])
genes = df["Name"].values
mat = floor_tpm(df[tissues].values.astype(float))
n_tissues, n_genes = len(tissues), mat.shape[0]

rng = np.random.default_rng(SEED)

if RESAMPLE_MODE == "cluster":
    Z = pd.read_csv(IN_LINKAGE).values
    labels = fcluster(Z, t=N_BLOCKS, criterion="maxclust")
    blocks = [np.where(labels == b)[0] for b in np.unique(labels)]
    for b, idx in zip(np.unique(labels), blocks):
        print(f"block {b} ({len(idx)}): {[tissues[i] for i in idx]}", file=sys.stderr)

reps = np.full((N_REPS, n_genes), np.nan)
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
    reps[r] = np.where(np.isfinite(tau_r), tau_r, np.nan)

    if (r + 1) % max(1, N_REPS // 10) == 0:
        print(f"{TAU_TYPE}/{RESAMPLE_MODE}: rep {r + 1}/{N_REPS}", file=sys.stderr)

out = pd.DataFrame({
    "gene": genes,
    f"tau_{TAU_TYPE}_{RESAMPLE_MODE}_var": np.nanvar(reps, axis=0, ddof=1),
    f"tau_{TAU_TYPE}_{RESAMPLE_MODE}_mean": np.nanmean(reps, axis=0),
})
out.to_csv(OUT_CSV, index=False)
print(f"genes={n_genes} reps={N_REPS} tau_type={TAU_TYPE} mode={RESAMPLE_MODE}", file=sys.stderr)
