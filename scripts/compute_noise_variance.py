"""Monte Carlo noise-robustness test: perturb log10(TPM+1) with iid N(0, SIGMA) noise
on every currently-nonzero (TPM>=1, post floor_tpm) entry of the full base_gct --
entries already floored to 0 (TPM<1) are left untouched, and noisy values that drop
below the TPM>=1 threshold are NOT re-floored (see gls_tau.apply_nonzero_log_noise).
The weight vector (IN_WEIGHTS_CSV) is frozen at its already-fitted value -- not refit
per replicate, since it's an aggregate statistic over the whole gene set and
insensitive to per-cell noise at this scale (see applying_noise.md). Reports the plain
Monte Carlo variance per gene across replicates (np.var, no jackknife (n-1)/n
correction -- these are iid noise draws, not delete-1 subsamples)."""
import sys
import numpy as np
import pandas as pd
from gls_tau import floor_tpm, apply_nonzero_log_noise, weighted_tau

FULL_GCT = sys.argv[1]
IN_WEIGHTS_CSV = sys.argv[2]
WEIGHT_COL = sys.argv[3]
METHOD_NAME = sys.argv[4]
N_REPLICATES = int(sys.argv[5])
SIGMA = float(sys.argv[6])
SEED = int(sys.argv[7])
OUT_SUMMARY_CSV = sys.argv[8]
OUT_LONG_CSV = sys.argv[9]

df = pd.read_csv(FULL_GCT, sep="\t", skiprows=2)
tissues = list(df.columns[2:])
genes = df["Name"].values
n_genes = len(genes)
mat = floor_tpm(df[tissues].values.astype(float))  # TPM<1 -> 0, clean baseline

weights_df = pd.read_csv(IN_WEIGHTS_CSV).set_index("tissue")
w = weights_df.loc[tissues, WEIGHT_COL].values

rng = np.random.default_rng(SEED)
tau_u_mat = np.full((n_genes, N_REPLICATES), np.nan)
tau_w_mat = np.full((n_genes, N_REPLICATES), np.nan)

for rep_i in range(N_REPLICATES):
    mat_noisy = apply_nonzero_log_noise(mat, SIGMA, rng)
    tau_u_mat[:, rep_i] = weighted_tau(mat_noisy, np.ones(len(tissues)))
    tau_w_mat[:, rep_i] = weighted_tau(mat_noisy, w)

summary = pd.DataFrame({
    "gene": genes,
    "tau_unweighted_noise_mean": np.nanmean(tau_u_mat, axis=1),
    "tau_unweighted_noise_var": np.nanvar(tau_u_mat, axis=1),
    f"tau_{METHOD_NAME}_weighted_noise_mean": np.nanmean(tau_w_mat, axis=1),
    f"tau_{METHOD_NAME}_weighted_noise_var": np.nanvar(tau_w_mat, axis=1),
})
summary.to_csv(OUT_SUMMARY_CSV, index=False)

long_df = pd.DataFrame({
    "gene": np.tile(genes, N_REPLICATES),
    "replicate": np.repeat(np.arange(N_REPLICATES), n_genes),
    "tau_unweighted_noisy": tau_u_mat.T.flatten(),
    f"tau_{METHOD_NAME}_weighted_noisy": tau_w_mat.T.flatten(),
})
long_df.to_csv(OUT_LONG_CSV, index=False)

print(f"genes={n_genes} tissues={len(tissues)} replicates={N_REPLICATES} sigma={SIGMA} "
      f"median_var_unweighted={np.nanmedian(summary['tau_unweighted_noise_var']):.6g} "
      f"median_var_{METHOD_NAME}={np.nanmedian(summary[f'tau_{METHOD_NAME}_weighted_noise_var']):.6g}",
      file=sys.stderr)
