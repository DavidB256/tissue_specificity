"""Bootstrap variance of tau_unweighted / tau_weighted_gls_ridge per gene, under two
resampling schemes:
  tissue  -- resample the 68 GTEx tissues with replacement (standard tissue bootstrap).
  cluster -- resample 11 top-level correlation-tree tissue blocks with replacement
             (cut at t=11, i.e. first 10 splits of the dendrogram); each drawn block
             contributes all its member tissues. 9 blocks are hardcoded per the given
             spec (testis; EBV-lymphocytes+whole blood; liver x3; Brain_*; Pancreas_*;
             muscle/heart x3; fibroblasts; SI-lymphoid+spleen; pituitary); the
             remaining tissues (not in any of those) are split into the 2 remaining
             large blocks by average-linkage clustering on 1-corr distance, matching
             what an actual t=11 cut of the full tree gives for this leftover set.
             Tests whether GLS/ridge weighting is specifically robust to whole
             correlated blocks (e.g. all 13 Brain_* tissues) being duplicated/dropped
             together, vs single tissues under plain tissue resampling.

Seed is fixed (not a CLI arg) so tissue/cluster draws are identical across the
unweighted and ridge jobs for the same mode -- keeps the two tau types directly
comparable per rep even though they run as separate parallel snakemake jobs."""
import sys
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from gls_tau import tissue_corr, ridge_lambda_search, weighted_tau, floor_tpm

IN_GCT = sys.argv[1]
IN_CORR = sys.argv[2]
TAU_TYPE = sys.argv[3]       # unweighted | ridge
RESAMPLE_MODE = sys.argv[4]  # tissue | cluster
N_REPS = int(sys.argv[5])
OUT_CSV = sys.argv[6]

SEED = 0

EXPLICIT_BLOCKS = [
    ["Testis"],
    ["Cells_EBV-transformed_lymphocytes", "Whole_Blood"],
    ["Liver", "Liver_Hepatocyte", "Liver_Mixed_Cell"],
    None,  # Brain_* -- filled in below
    None,  # Pancreas* -- filled in below
    ["Muscle_Skeletal", "Heart_Atrial_Appendage", "Heart_Left_Ventricle"],
    ["Cells_Cultured_fibroblasts"],
    ["Small_Intestine_Terminal_Ileum_Lymphoid_Aggregate", "Spleen"],
    ["Pituitary"],
]


def build_blocks(tissues, corr_df):
    blocks = [b for b in EXPLICIT_BLOCKS if b is not None]
    blocks.append([t for t in tissues if t.startswith("Brain_")])
    blocks.append([t for t in tissues if t.startswith("Pancreas")])
    named = set(t for b in blocks for t in b)
    remaining = [t for t in tissues if t not in named]
    sub = corr_df.loc[remaining, remaining].values
    dist = 1 - sub
    np.fill_diagonal(dist, 0.0)
    dist = (dist + dist.T) / 2
    Z = linkage(squareform(dist, checks=False), method="average")
    labels = fcluster(Z, t=2, criterion="maxclust")
    for b in np.unique(labels):
        blocks.append([remaining[i] for i in np.where(labels == b)[0]])
    tissue_to_idx = {t: i for i, t in enumerate(tissues)}
    idx_blocks = [np.array([tissue_to_idx[t] for t in b]) for b in blocks]
    for b in blocks:
        print(f"block ({len(b)}): {b}", file=sys.stderr)
    return idx_blocks


df = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(df.columns[2:])
genes = df["Name"].values
mat = floor_tpm(df[tissues].values.astype(float))
n_tissues, n_genes = len(tissues), mat.shape[0]

rng = np.random.default_rng(SEED)

if RESAMPLE_MODE == "cluster":
    corr_df = pd.read_csv(IN_CORR, index_col=0)
    blocks = build_blocks(tissues, corr_df)
    print(f"n_blocks={len(blocks)}", file=sys.stderr)

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
