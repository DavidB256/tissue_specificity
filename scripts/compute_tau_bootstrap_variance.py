"""Bootstrap variance of tau_unweighted / tau_weighted_gls_ridge per gene, under three
resampling schemes:
  tissue     -- resample the 68 GTEx tissues with replacement (standard tissue
                bootstrap).
  cluster10  -- resample the 10 tissue blocks (t=10 dendrogram cut, see
                identify_tissue_clusters.py) with replacement; each drawn block
                contributes all its member tissues.
  cluster20  -- same, with a finer t=20 dendrogram cut (20 tissue blocks).
Cluster resampling tests whether GLS/ridge weighting is specifically robust to whole
correlated blocks (e.g. all 13 Brain_* tissues) being duplicated/dropped together, vs
single tissues under plain tissue resampling.

Seed is fixed (not a CLI arg) so tissue/cluster draws are identical across the
unweighted and ridge jobs for the same mode -- keeps the two tau types directly
comparable per rep even though they run as separate parallel snakemake jobs."""
import sys
import numpy as np
import pandas as pd
from gls_tau import tissue_corr, tissue_corr_expressed_only, ridge_lambda_search, weighted_tau, floor_tpm

IN_GCT = sys.argv[1]
IN_CLUSTERS = sys.argv[2]
TAU_TYPE = sys.argv[3]       # unweighted | ridge
RESAMPLE_MODE = sys.argv[4]  # tissue | cluster10 | cluster20
N_REPS = int(sys.argv[5])
OUT_CSV = sys.argv[6]
CORR_MODE = sys.argv[7] if len(sys.argv) > 7 else "all"  # all | expressed (TPM>=1 in >=1 tissue)
corr_fn = tissue_corr_expressed_only if CORR_MODE == "expressed" else tissue_corr

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
