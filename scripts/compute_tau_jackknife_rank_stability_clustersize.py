"""Jackknife rank-stability for the cluster-size-downweighted tau variant (see
compute_tau_clustersize.py). Dropped tissue(s) determined by RESAMPLE_MODE (same
tissue/cluster10/cluster20 schemes as compute_tau_jackknife_rank_stability.py); tau
weighting always uses the fixed k-block assignment from tissue_clusters_{k}.csv --
block LABELS are static, block SIZE is recomputed per rep from tissues actually
present after the drop."""
import sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from gls_tau import weighted_tau, floor_tpm, gene_subset_mask

IN_GCT = sys.argv[1]
IN_TAU = sys.argv[2]                 # gene_specificity_clustersize{k}/weighted_tau.csv
IN_LEAVEOUT_CLUSTERS = sys.argv[3]   # leave-out block definition (RESAMPLE_MODE-dependent)
IN_WEIGHT_CLUSTERS = sys.argv[4]     # fixed k-block assignment used for weighting
RESAMPLE_MODE = sys.argv[5]          # tissue | cluster10 | cluster20
OUT_CSV = sys.argv[6]
GENE_SUBSET = sys.argv[7] if len(sys.argv) > 7 else "all"

TAU_COL = "tau_clustersize"

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

weight_blocks = pd.read_csv(IN_WEIGHT_CLUSTERS).set_index("tissue").loc[tissues, "block"].values

if RESAMPLE_MODE == "tissue":
    units = tissues
    leave_out_idx = [np.array([i]) for i in range(n_tissues)]
else:
    clusters = pd.read_csv(IN_LEAVEOUT_CLUSTERS)
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
    blocks_r = weight_blocks[idx]
    uniq, counts = np.unique(blocks_r, return_counts=True)
    size_map = dict(zip(uniq, counts))
    w = np.array([1.0 / size_map[b] for b in blocks_r])
    tau_r = weighted_tau(mat_r, w)

    ok = np.isfinite(tau_ref) & np.isfinite(tau_r) & subset_mask
    rho, _ = spearmanr(tau_ref[ok], tau_r[ok])
    rhos[r] = rho
    n_used[r] = ok.sum()
    print(f"clustersize/{RESAMPLE_MODE}: rep {r + 1}/{n_reps} (drop {units[r]})", file=sys.stderr)

out = pd.DataFrame({
    "rep": np.arange(n_reps),
    "unit": units,
    f"spearman_rho_clustersize_{RESAMPLE_MODE}": rhos,
    "n_genes": n_used,
})
out.to_csv(OUT_CSV, index=False)
print(f"reps={n_reps} mode={RESAMPLE_MODE} median_rho={np.nanmedian(rhos):.4f}", file=sys.stderr)
