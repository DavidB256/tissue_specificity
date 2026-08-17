"""Stability tests for the tau variants:
(1) pairwise Spearman rank correlation across the whole gene set.
(2) leave-one-tissue-out (LOTO): recompute all three tau variants with each of the
    68 tissues dropped in turn, then report each gene's std across the 68 reruns.
    Weighted variants should show *lower* LOTO std if the GLS correction is actually
    discounting genuinely redundant tissues -- dropping one member of a correlated
    cluster (e.g. one brain subregion) should barely move a gene's specificity score
    once that redundancy is accounted for."""
import sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from gls_tau import tissue_corr_expressed_only, gls_weights, ridge_lambda_search, weighted_tau, floor_tpm, censor_floor

IN_GCT = sys.argv[1]
IN_TAU = sys.argv[2]
OUT_SPEARMAN_CSV = sys.argv[3]
OUT_LOTO_CSV = sys.argv[4]

VARIANTS = ["tau_unweighted", "tau_weighted_gls", "tau_weighted_gls_ridge"]

# --- (1) Spearman rank correlation between variants, whole gene set ---
tau_df = pd.read_csv(IN_TAU)[["gene"] + VARIANTS].dropna()
corr = pd.DataFrame(index=VARIANTS, columns=VARIANTS, dtype=float)
for a in VARIANTS:
    for b in VARIANTS:
        corr.loc[a, b] = spearmanr(tau_df[a], tau_df[b]).correlation
corr.to_csv(OUT_SPEARMAN_CSV)
print(f"Spearman correlation between variants (n={len(tau_df)} genes):\n{corr}", file=sys.stderr)

# --- (2) leave-one-tissue-out reproducibility ---
df = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(df.columns[2:])
genes = df["Name"].values
mat_full = floor_tpm(df[tissues].values.astype(float))  # TPM<1 -> 0
n_tissues, n_genes = len(tissues), mat_full.shape[0]

loto = {v: np.full((n_tissues, n_genes), np.nan) for v in VARIANTS}

for k in range(n_tissues):
    mat = mat_full[:, np.arange(n_tissues) != k]
    n = mat.shape[1]
    R = tissue_corr_expressed_only(mat)
    w_unreg = gls_weights(R, 0.0)
    lam_star, w_reg = ridge_lambda_search(R)

    loto["tau_unweighted"][k] = weighted_tau(mat, np.ones(n))
    loto["tau_weighted_gls"][k] = censor_floor(weighted_tau(mat, w_unreg))
    loto["tau_weighted_gls_ridge"][k] = weighted_tau(mat, w_reg)
    print(f"LOTO {k + 1}/{n_tissues} ({tissues[k]}) done, lambda*={lam_star:.4g}", file=sys.stderr)

out = pd.DataFrame({"gene": genes})
for v in VARIANTS:
    out[f"{v}_loto_std"] = np.nanstd(loto[v], axis=0)
    out[f"{v}_loto_mean"] = np.nanmean(loto[v], axis=0)
out.to_csv(OUT_LOTO_CSV, index=False)

summary = out[[f"{v}_loto_std" for v in VARIANTS]].median()
print(f"median per-gene LOTO std across variants:\n{summary}", file=sys.stderr)
