"""Parametric bootstrap for chess_tx_protein_coding, standing in for
bootstrap_gtex_v11_tx_protein_coding.py's real donor resampling -- the CHESS
bodysite raw table has no per-donor values, only per-(transcript, tissue)
summary stats (min/25/50/75/max/count, from reformat_chess_tx_tissue_stats.py).

Per resample, per (transcript, tissue) cell: treat the 5 reported quantiles as
a piecewise-linear approximation of that cell's quantile function Q(u) (u in
[0,1]), and simulate "the median of `count` iid draws from Q" directly via the
order-statistic shortcut -- for n iid draws of a continuous distribution, the
median equals Q(U_(k)) where U_(k) is the k-th order statistic of n iid
Uniform(0,1) draws (k = ceil(n/2)), and U_(k) ~ Beta(k, n-k+1) exactly. So one
Beta draw per cell per resample, pushed through the piecewise-linear Q, gives
the resampled per-tissue value -- no need to actually draw `count` individual
values. A cell with count==0 (tissue never observed for that transcript) is
deterministic 0, matching build_chess_tx_gct.py's convention.

Caveat vs. the real bootstrap: cells are resampled INDEPENDENTLY across
transcripts and tissues (we only have marginal per-cell quantiles, no
donor-level joint data), whereas the real script resamples the same donor
index set across all transcripts in a tissue, preserving cross-transcript
covariance within a donor. Gene-level TPM here is therefore a sum of
independently-resampled transcripts, not a true per-donor co-resample --
likely UNDERSTATES gene-level bootstrap variance relative to the real method.

Universe, outputs, and tau computation (tau_unweighted, tau_weighted with the
already-fitted chess_tx_protein_coding free_cluster_weights_max/hierclust0.1
weight vector) are otherwise identical in shape to
bootstrap_gtex_v11_tx_protein_coding.py -- see that script's docstring."""
import sys
import numpy as np
import pandas as pd
from scipy import sparse
from gls_tau import weighted_tau

IN_PC_GCT = sys.argv[1]
IN_TISSUE_STATS_CSV = sys.argv[2]
IN_WEIGHTS_CSV = sys.argv[3]
OUT_DIR = sys.argv[4]
N_BOOTSTRAP = int(sys.argv[5]) if len(sys.argv) > 5 else 200
SEED = int(sys.argv[6]) if len(sys.argv) > 6 else 0

STATS = ["min", "25", "50", "75", "max", "count"]
PROBS = np.array([0.0, 0.25, 0.5, 0.75, 1.0])

rng = np.random.default_rng(SEED)

# ---- 1. "expressed" universe: same tx_total>0 filter as compute_tx_info_content.py ----
pc = pd.read_csv(IN_PC_GCT, sep="\t", skiprows=2)
tissues_pc = list(pc.columns[2:])
mat_pc = pc[tissues_pc].to_numpy(dtype=float)
tx_total_full = mat_pc.sum(axis=1)
keep = tx_total_full > 0
tx_ids = pc["Name"].to_numpy()[keep]
tx_gene_ids = pc["Description"].to_numpy()[keep]
n_tx = len(tx_ids)
gene_ids = np.array(sorted(pd.unique(tx_gene_ids)))
n_gene = len(gene_ids)
gene_index = {g: i for i, g in enumerate(gene_ids)}
gene_row_of_tx = np.array([gene_index[g] for g in tx_gene_ids])
G = sparse.csr_matrix((np.ones(n_tx), (gene_row_of_tx, np.arange(n_tx))), shape=(n_gene, n_tx))
print(f"universe: {n_tx} expressed transcripts, {n_gene} expressed genes "
      f"(of {len(pc)} total transcripts)", file=sys.stderr)
del pc, mat_pc, tx_total_full, keep

# ---- 2. per-cell quantiles/count, reindexed to our expressed tx x tissues_pc order ----
weights_df = pd.read_csv(IN_WEIGHTS_CSV).set_index("tissue")
tissues = tissues_pc
assert set(tissues) == set(weights_df.index), "tissue set mismatch vs weights.csv"
w = weights_df.loc[tissues, "weight"].to_numpy()
n_tissues = len(tissues)

stats_cols = [f"{t}_{s}" for t in tissues for s in STATS]
stats = pd.read_csv(IN_TISSUE_STATS_CSV, usecols=["transcript_id"] + stats_cols) \
    .set_index("transcript_id").reindex(tx_ids)

q0 = stats[[f"{t}_min" for t in tissues]].to_numpy(dtype=np.float64)
q25 = stats[[f"{t}_25" for t in tissues]].to_numpy(dtype=np.float64)
q50 = stats[[f"{t}_50" for t in tissues]].to_numpy(dtype=np.float64)
q75 = stats[[f"{t}_75" for t in tissues]].to_numpy(dtype=np.float64)
q100 = stats[[f"{t}_max" for t in tissues]].to_numpy(dtype=np.float64)
count = stats[[f"{t}_count" for t in tissues]].to_numpy(dtype=np.float64)
count = np.nan_to_num(count, nan=0.0)
observed = count > 0
print(f"observed (transcript, tissue) cells: {observed.sum()} / {observed.size}", file=sys.stderr)
del stats


def piecewise_quantile(u, vals):
    """vals: list of 5 (n_tx, n_tissues) arrays at PROBS=[0,.25,.5,.75,1]."""
    stacked = np.stack(vals, axis=-1)  # (n_tx, n_tissues, 5)
    seg = np.clip(np.searchsorted(PROBS, u, side="right") - 1, 0, 3)
    lo_p, hi_p = PROBS[seg], PROBS[seg + 1]
    lo_v = np.take_along_axis(stacked, seg[..., None], axis=-1)[..., 0]
    hi_v = np.take_along_axis(stacked, (seg + 1)[..., None], axis=-1)[..., 0]
    frac = (u - lo_p) / (hi_p - lo_p)
    return lo_v + frac * (hi_v - lo_v)


# Beta(k, n-k+1) shape params for the median order statistic (k = ceil(n/2)); placeholder
# (1, 1) on unobserved cells avoids invalid-parameter draws, overwritten with 0 below.
k = np.where(observed, np.ceil(count / 2.0), 1.0)
beta_a = k
beta_b = np.where(observed, count - k + 1.0, 1.0)

# ---- 3. Welford accumulators ----
def welford_init(shape):
    return np.zeros(shape), np.zeros(shape), np.zeros(shape, dtype=int)


tx_tpm_mean, tx_tpm_M2, tx_tpm_n = welford_init((n_tx, n_tissues))
gene_tpm_mean, gene_tpm_M2, gene_tpm_n = welford_init((n_gene, n_tissues))

tx_tauU_mean, tx_tauU_M2, tx_tauU_n = welford_init(n_tx)
tx_tauW_mean, tx_tauW_M2, tx_tauW_n = welford_init(n_tx)
gene_tauU_mean, gene_tauU_M2, gene_tauU_n = welford_init(n_gene)
gene_tauW_mean, gene_tauW_M2, gene_tauW_n = welford_init(n_gene)

count_valid_U = np.zeros(n_tx, dtype=int)
count_gt_U = np.zeros(n_tx, dtype=int)
count_valid_W = np.zeros(n_tx, dtype=int)
count_gt_W = np.zeros(n_tx, dtype=int)


def welford_update(mean, M2, n, x, mask=None):
    if mask is None:
        n += 1
        delta = x - mean
        mean += delta / n
        M2 += delta * (x - mean)
    else:
        n[mask] += 1
        delta = np.where(mask, x - mean, 0.0)
        mean[mask] += delta[mask] / n[mask]
        delta2 = np.where(mask, x - mean, 0.0)
        M2[mask] += (delta * delta2)[mask]


# ---- 4. bootstrap loop ----
for b in range(N_BOOTSTRAP):
    u = rng.beta(beta_a, beta_b)
    tx_mat_b = piecewise_quantile(u, [q0, q25, q50, q75, q100])
    tx_mat_b = np.where(observed, tx_mat_b, 0.0)

    gene_mat_b = np.asarray(G @ tx_mat_b)

    tau_txU = weighted_tau(tx_mat_b, np.ones(n_tissues))
    tau_txW = weighted_tau(tx_mat_b, w)
    tau_geU = weighted_tau(gene_mat_b, np.ones(n_tissues))
    tau_geW = weighted_tau(gene_mat_b, w)

    welford_update(tx_tpm_mean, tx_tpm_M2, tx_tpm_n, tx_mat_b)
    welford_update(gene_tpm_mean, gene_tpm_M2, gene_tpm_n, gene_mat_b)

    mU = np.isfinite(tau_txU)
    welford_update(tx_tauU_mean, tx_tauU_M2, tx_tauU_n, tau_txU, mU)
    mW = np.isfinite(tau_txW)
    welford_update(tx_tauW_mean, tx_tauW_M2, tx_tauW_n, tau_txW, mW)
    mgU = np.isfinite(tau_geU)
    welford_update(gene_tauU_mean, gene_tauU_M2, gene_tauU_n, tau_geU, mgU)
    mgW = np.isfinite(tau_geW)
    welford_update(gene_tauW_mean, gene_tauW_M2, gene_tauW_n, tau_geW, mgW)

    tau_geU_aligned = tau_geU[gene_row_of_tx]
    validU = mU & mgU[gene_row_of_tx]
    count_valid_U += validU
    count_gt_U += validU & (tau_txU > tau_geU_aligned)

    tau_geW_aligned = tau_geW[gene_row_of_tx]
    validW = mW & mgW[gene_row_of_tx]
    count_valid_W += validW
    count_gt_W += validW & (tau_txW > tau_geW_aligned)

    if (b + 1) % 10 == 0:
        print(f"  bootstrap {b + 1}/{N_BOOTSTRAP} done", file=sys.stderr)

# ---- 5. finalize + write ----
def se(M2, n):
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(n > 1, np.sqrt(M2 / (n - 1)), np.nan)


tx_tpm_se = se(tx_tpm_M2, tx_tpm_n)
gene_tpm_se = se(gene_tpm_M2, gene_tpm_n)

pd.DataFrame(tx_tpm_se, columns=tissues).assign(transcript_id=tx_ids, gene_id=tx_gene_ids) \
    [["transcript_id", "gene_id"] + tissues].to_csv(f"{OUT_DIR}/tx_tpm_se.csv", index=False)
pd.DataFrame(gene_tpm_se, columns=tissues).assign(gene_id=gene_ids)[["gene_id"] + tissues] \
    .to_csv(f"{OUT_DIR}/gene_tpm_se.csv", index=False)

tx_summary = pd.DataFrame({
    "transcript_id": tx_ids,
    "gene_id": tx_gene_ids,
    "se_tau_unweighted": se(tx_tauU_M2, tx_tauU_n),
    "se_tau_weighted": se(tx_tauW_M2, tx_tauW_n),
    "prob_tx_tau_gt_gene_unweighted": np.where(count_valid_U > 0, count_gt_U / np.maximum(count_valid_U, 1), np.nan),
    "prob_tx_tau_gt_gene_weighted": np.where(count_valid_W > 0, count_gt_W / np.maximum(count_valid_W, 1), np.nan),
    "n_valid_unweighted": count_valid_U,
    "n_valid_weighted": count_valid_W,
})
tx_summary.to_csv(f"{OUT_DIR}/tx_tau_summary.csv", index=False)

gene_summary = pd.DataFrame({
    "gene_id": gene_ids,
    "se_tau_unweighted": se(gene_tauU_M2, gene_tauU_n),
    "se_tau_weighted": se(gene_tauW_M2, gene_tauW_n),
})
gene_summary.to_csv(f"{OUT_DIR}/gene_tau_summary.csv", index=False)

print(f"N_BOOTSTRAP={N_BOOTSTRAP} seed={SEED} -> {OUT_DIR}/{{tx,gene}}_tpm_se.csv, "
      f"{{tx,gene}}_tau_summary.csv", file=sys.stderr)
