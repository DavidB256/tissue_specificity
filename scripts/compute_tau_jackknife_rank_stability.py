"""Jackknife counterpart to compute_tau_bootstrap_rank_stability.py: instead of
resampling-with-replacement, deterministically leave one unit out at a time and
measure Spearman rho between the leave-one-out tau ranking and the canonical
(full-data) tau ranking from weighted_tau.csv.
  tissue     -- leave one of the 68 GTEx tissues out (68 reps).
  cluster10  -- leave one whole t=10 tissue block out (10 reps).
  cluster20  -- leave one whole t=20 tissue block out (20 reps).
No seed/N_REPS needed: jackknife is exhaustive over leave-out units."""
import sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from gls_tau import tissue_corr, ridge_lambda_search, weighted_tau, floor_tpm, gene_subset_mask

IN_GCT = sys.argv[1]
IN_TAU = sys.argv[2]
IN_CLUSTERS = sys.argv[3]
TAU_TYPE = sys.argv[4]       # unweighted | ridge
RESAMPLE_MODE = sys.argv[5]  # tissue | cluster10 | cluster20
OUT_CSV = sys.argv[6]
GENE_SUBSET = sys.argv[7] if len(sys.argv) > 7 else "all"  # all | highspec | lowspec | multitissue | onlyneural | maxneural

TAU_COL = "tau_unweighted" if TAU_TYPE == "unweighted" else "tau_weighted_gls_ridge"

df = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(df.columns[2:])
genes = df["Name"].values
mat = floor_tpm(df[tissues].values.astype(float))

expressed = mat.max(axis=1) > 0
print(f"dropping {(~expressed).sum()}/{len(genes)} genes with no tissue TPM>=1", file=sys.stderr)
genes, mat = genes[expressed], mat[expressed]
n_tissues = len(tissues)

tau_df = pd.read_csv(IN_TAU).set_index("gene").loc[genes]
tau_ref = tau_df[TAU_COL].values
subset_mask = gene_subset_mask(GENE_SUBSET, mat, tissues, tau_df["tau_unweighted"].values)
print(f"gene_subset={GENE_SUBSET}: {subset_mask.sum()}/{len(genes)} genes", file=sys.stderr)

if RESAMPLE_MODE == "tissue":
    units = tissues
    leave_out_idx = [np.array([i]) for i in range(n_tissues)]
else:
    clusters = pd.read_csv(IN_CLUSTERS)
    tissue_to_idx = {t: i for i, t in enumerate(tissues)}
    blocks = [np.array([tissue_to_idx[t] for t in g["tissue"]])
              for _, g in clusters.groupby("block")]
    units = [f"block_{b}" for b, _ in clusters.groupby("block")]
    leave_out_idx = blocks

n_reps = len(units)
rhos = np.full(n_reps, np.nan)
n_used = np.full(n_reps, 0)
for r, drop_idx in enumerate(leave_out_idx):
    idx = np.setdiff1d(np.arange(n_tissues), drop_idx)
    mat_r = mat[:, idx]
    if TAU_TYPE == "unweighted":
        tau_r = weighted_tau(mat_r, np.ones(len(idx)))
    else:
        R = tissue_corr(mat_r)
        _, w = ridge_lambda_search(R)
        tau_r = weighted_tau(mat_r, w)

    ok = np.isfinite(tau_ref) & np.isfinite(tau_r) & subset_mask
    rho, _ = spearmanr(tau_ref[ok], tau_r[ok])
    rhos[r] = rho
    n_used[r] = ok.sum()
    print(f"{TAU_TYPE}/{RESAMPLE_MODE}: rep {r + 1}/{n_reps} (drop {units[r]})", file=sys.stderr)

out = pd.DataFrame({
    "rep": np.arange(n_reps),
    "unit": units,
    f"spearman_rho_{TAU_TYPE}_{RESAMPLE_MODE}": rhos,
    "n_genes": n_used,
})
out.to_csv(OUT_CSV, index=False)
print(f"reps={n_reps} tau_type={TAU_TYPE} mode={RESAMPLE_MODE} "
      f"median_rho={np.nanmedian(rhos):.4f}", file=sys.stderr)
