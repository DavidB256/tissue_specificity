"""Nonparametric multi-sample bootstrap for gtex_v11_tx_protein_coding: resample
GTEx donor samples independently within each tissue (with replacement, same n
per tissue as the real data), recompute per-tissue TPM (median across the
resampled donors, NOT floored -- unlike the real pipeline, so a resample
doesn't lose a lowly-expressed tx to the TPM<1 floor), then gene-level TPM
(sum of its kept transcripts' resampled TPM per tissue), then tau_unweighted
(all-ones tissue weights) and tau_weighted (the ALREADY-FITTED shared
free_cluster_weights hierclust0.1/max-agg vector, IN_WEIGHTS_CSV -- frozen,
not refit per resample) for both resolutions. Repeated N_BOOTSTRAP times.

Universe: restricted to the "expressed" transcripts/genes already established
in compute_tx_info_content.py (nonzero total TPM in the real, floored,
median-aggregated data) -- an all-zero-in-the-real-data tx/gene is excluded
outright, never considered across resamples.

Outputs (OUT_DIR):
  tx_tpm_se.csv / gene_tpm_se.csv -- bootstrap SE of per-tissue TPM, one row
    per tx/gene, one column per tissue (plus gene_id for tx file).
  tx_tau_summary.csv -- transcript_id, gene_id, se_tau_unweighted,
    se_tau_weighted, prob_tx_tau_gt_gene_unweighted,
    prob_tx_tau_gt_gene_weighted, n_valid_unweighted, n_valid_weighted (the
    last two: resamples where both the tx's and its gene's tau were defined,
    i.e. not all-zero in that resample -- the denominator for the prob cols).
  gene_tau_summary.csv -- gene_id, se_tau_unweighted, se_tau_weighted.

Welford's online algorithm throughout (mean/M2 running accumulators) so the
B resampled matrices are never all held in memory at once -- only the
(n_tx or n_gene) x n_tissues-shaped running accumulators."""
import sys
import numpy as np
import pandas as pd
from scipy import sparse
from gls_tau import weighted_tau

IN_RAW_TPM = sys.argv[1]
IN_SAMPLE_TISSUE_CSV = sys.argv[2]
IN_PC_GCT = sys.argv[3]
IN_WEIGHTS_CSV = sys.argv[4]
OUT_DIR = sys.argv[5]
N_BOOTSTRAP = int(sys.argv[6]) if len(sys.argv) > 6 else 200
SEED = int(sys.argv[7]) if len(sys.argv) > 7 else 0

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

# ---- 2. sample -> tissue map, tissue order matching weights.csv ----
sample_tissue = pd.read_csv(IN_SAMPLE_TISSUE_CSV).set_index("sample_id")["tissue"]
weights_df = pd.read_csv(IN_WEIGHTS_CSV).set_index("tissue")

header = pd.read_csv(IN_RAW_TPM, sep="\t", nrows=0).columns
sample_cols = [c for c in header[2:] if c in sample_tissue.index]
tissues = sorted(sample_tissue.loc[sample_cols].unique())
assert set(tissues) == set(weights_df.index), "tissue set mismatch vs weights.csv"
w = weights_df.loc[tissues, "weight"].to_numpy()
n_tissues = len(tissues)
sample_col_pos = {c: i for i, c in enumerate(sample_cols)}
tissue_positions = {t: np.array([sample_col_pos[c] for c in sample_cols if sample_tissue[c] == t])
                     for t in tissues}
print(f"tissues={n_tissues}, sample columns={len(sample_cols)}", file=sys.stderr)

# ---- 3. stream the raw per-sample TPM table once, keep only our expressed tx rows ----
tx_id_set = set(tx_ids)
CHUNKSIZE = 4000
dtype = {"transcript_id": str, "gene_id": str}
dtype.update({c: np.float32 for c in sample_cols})
collected_rows = []
collected_ids = []
n_rows_scanned = 0
reader = pd.read_csv(IN_RAW_TPM, sep="\t", usecols=["transcript_id"] + sample_cols,
                      dtype=dtype, chunksize=CHUNKSIZE)
for chunk in reader:
    hit = chunk[chunk["transcript_id"].isin(tx_id_set)]
    if len(hit):
        collected_rows.append(hit[sample_cols].to_numpy(dtype=np.float32))
        collected_ids.append(hit["transcript_id"].to_numpy())
    n_rows_scanned += len(chunk)
    if n_rows_scanned % 40000 == 0:
        print(f"  ...{n_rows_scanned} raw rows scanned, {sum(len(x) for x in collected_ids)} matched",
              file=sys.stderr)

raw_ids = np.concatenate(collected_ids)
raw_mat = np.concatenate(collected_rows, axis=0)
del collected_rows, collected_ids
print(f"matched {len(raw_ids)} / {n_tx} expressed transcripts in raw table", file=sys.stderr)

# reorder raw_mat rows to match our canonical tx_ids order
pos = pd.Series(np.arange(len(raw_ids)), index=raw_ids)
reorder = pos.loc[tx_ids].to_numpy()
tx_raw_mat = raw_mat[reorder]
del raw_mat, raw_ids, pos, reorder
print(f"tx_raw_mat shape={tx_raw_mat.shape} ({tx_raw_mat.nbytes / 1e9:.2f} GB)", file=sys.stderr)

# ---- 4. Welford accumulators ----
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


# ---- 5. bootstrap loop ----
for b in range(N_BOOTSTRAP):
    tx_mat_b = np.empty((n_tx, n_tissues), dtype=np.float64)
    for ti, t in enumerate(tissues):
        idx = tissue_positions[t]
        resampled = rng.choice(idx, size=len(idx), replace=True)
        tx_mat_b[:, ti] = np.median(tx_raw_mat[:, resampled], axis=1)

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

# ---- 6. finalize + write ----
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
