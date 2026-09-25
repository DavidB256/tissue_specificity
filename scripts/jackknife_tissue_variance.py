"""Delete-one-TISSUE jackknife variance of tau_unweighted and tau_free_cluster_weighted,
parallel to jackknife_cluster_variance.py but ignoring tissue correlation structure: drops
one raw tissue at a time (not a whole hierarchical cluster) from the full-resolution panel.
free_cluster_weights fitting is unchanged (still one weight per hierarchical cluster,
CLUSTERS_CSV, DIST_THRESHOLD=0.1) -- only the jackknife deletion unit differs. A cluster's
weight dimension is dropped only if the deleted tissue was its last surviving member.
Per replicate:
  - tau_unweighted recomputed directly on the reduced (n-1 tissue) panel
  - the leg's collapsed reference is rebuilt from the reduced panel (same mean/max-of-MEMBERS
    logic as build_gct_onegroup.py) and free_cluster_weights is re-fit (differential_evolution)
    against it, warm-started at the already-discovered optimal weights (see
    jackknife_cluster_variance.py for the warm-start rationale)
Reports the standard delete-1 jackknife variance per gene: (n-1)/n * sum((theta_i-mean)^2).
Args: FULL_GCT CLUSTERS_CSV LEG AGG OPT_CLUSTER_WEIGHTS_CSV OUT_SUMMARY_CSV OUT_LONG_CSV
      MEMBERS... (or a single "prefix:XXX" arg, e.g. prefix:Brain, for prefix-matched legs)"""
import sys
import time
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from scipy.optimize import differential_evolution
from gls_tau import weighted_tau, floor_tpm

FULL_GCT = sys.argv[1]
CLUSTERS_CSV = sys.argv[2]
LEG = sys.argv[3]
AGG = sys.argv[4]
assert AGG in ("mean", "max"), AGG
OPT_CLUSTER_WEIGHTS_CSV = sys.argv[5]
OUT_SUMMARY_CSV = sys.argv[6]
OUT_LONG_CSV = sys.argv[7]
MEMBERS = sys.argv[8:]

SEED = 0
MAXITER = 300
POPSIZE = 15
JITTER_SD = 0.05

df = pd.read_csv(FULL_GCT, sep="\t", skiprows=2)
all_tissues = list(df.columns[2:])
genes = df["Name"].values
n_genes = len(genes)
raw_mat_full = df[all_tissues].values.astype(float)  # genes x tissues, unfloored

if len(MEMBERS) == 1 and MEMBERS[0].startswith("prefix:"):
    prefix = MEMBERS[0][len("prefix:"):]
    MEMBERS = [t for t in all_tissues if t.split("_")[0] == prefix]
missing = [m for m in MEMBERS if m not in all_tissues]
assert not missing, f"MEMBERS not in {FULL_GCT}: {missing}"

clusters = pd.read_csv(CLUSTERS_CSV).set_index("tissue").loc[all_tissues]
cluster_id_of_tissue = clusters["cluster_id"].values
n_tissues = len(all_tissues)
opt_weights = pd.read_csv(OPT_CLUSTER_WEIGHTS_CSV).set_index("cluster_id")["weight"]

print(f"leg={LEG} agg={AGG} tissues={n_tissues} clusters={clusters['cluster_id'].nunique()} "
      f"members={MEMBERS}", file=sys.stderr)

tau_u_mat = np.full((n_genes, n_tissues), np.nan)
tau_w_mat = np.full((n_genes, n_tissues), np.nan)
rho_log = []

for rep_i, dropped_tissue in enumerate(all_tissues):
    t0 = time.time()
    keep_mask = np.arange(n_tissues) != rep_i
    kept_tissues = [t for t, m in zip(all_tissues, keep_mask) if m]
    raw_mat_k = raw_mat_full[:, keep_mask]
    mat_k = floor_tpm(raw_mat_k)

    tau_u_k = weighted_tau(mat_k, np.ones(len(kept_tissues)))
    tau_u_mat[:, rep_i] = tau_u_k

    members_k = [m for m in MEMBERS if m in kept_tissues]
    other_k = [t for t in kept_tissues if t not in members_k]
    members_idx = [kept_tissues.index(m) for m in members_k]
    other_idx = [kept_tissues.index(t) for t in other_k]
    agg_fn = np.max if AGG == "max" else np.mean
    collapsed_col = agg_fn(raw_mat_k[:, members_idx], axis=1)
    ref_raw = np.concatenate([raw_mat_k[:, other_idx], collapsed_col[:, None]], axis=1)
    ref_mat = floor_tpm(ref_raw)
    tau_ref_k = weighted_tau(ref_mat, np.ones(ref_mat.shape[1]))

    remaining_clusters = sorted(set(cluster_id_of_tissue[keep_mask]))
    remap = {c: i for i, c in enumerate(remaining_clusters)}
    cluster_of_kept = np.array([remap[c] for c in cluster_id_of_tissue[keep_mask]])
    n_dims = len(remaining_clusters)

    valid = np.isfinite(tau_ref_k) & (mat_k.max(axis=1) > 0)
    mat_valid = mat_k[valid]
    tau_ref_valid = tau_ref_k[valid]
    max_i_valid = mat_valid.max(axis=1)
    argmax_idx_valid = mat_valid.argmax(axis=1)
    term_valid = 1 - mat_valid / max_i_valid[:, None]
    rank_ref = rankdata(tau_ref_valid)
    rank_ref_c = rank_ref - rank_ref.mean()
    denom_ref = np.sqrt((rank_ref_c ** 2).sum())

    def neg_spearman(x, cluster_of=cluster_of_kept, term=term_valid, argmax_idx=argmax_idx_valid,
                      rank_ref_c=rank_ref_c, denom_ref=denom_ref):
        w = x[cluster_of]
        tau = term @ w / (w.sum() - w[argmax_idx])
        rank_tau = rankdata(tau)
        rank_tau_c = rank_tau - rank_tau.mean()
        rho = (rank_tau_c @ rank_ref_c) / (np.sqrt((rank_tau_c ** 2).sum()) * denom_ref)
        return -rho

    x0 = np.array([opt_weights[c] for c in remaining_clusters])
    rng = np.random.default_rng(SEED)
    pop_size = POPSIZE * n_dims
    pop = np.clip(x0[None, :] + rng.normal(0, JITTER_SD, size=(pop_size, n_dims)), 0.01, 1.0)
    pop[0] = x0
    bounds = [(0.01, 1.0)] * n_dims

    result = differential_evolution(neg_spearman, bounds=bounds, seed=SEED, maxiter=MAXITER,
                                     tol=1e-8, polish=False, disp=False, init=pop)
    x_best = result.x
    best_rho = -result.fun
    w_best_kept = x_best[cluster_of_kept]
    tau_w_k = weighted_tau(mat_k, w_best_kept)
    tau_w_mat[:, rep_i] = tau_w_k

    rho_log.append({"dropped_tissue": dropped_tissue, "n_dims": n_dims, "best_rho": best_rho})
    print(f"[{rep_i + 1}/{n_tissues}] drop tissue={dropped_tissue} n_dims={n_dims} "
          f"best_rho={best_rho:.6f} ({time.time() - t0:.1f}s)", file=sys.stderr)


def jackknife_mean_var(mat):
    n_valid = np.sum(~np.isnan(mat), axis=1)
    mean = np.nanmean(mat, axis=1)
    with np.errstate(invalid="ignore"):
        sumsq = np.nansum((mat - mean[:, None]) ** 2, axis=1)
    var = np.where(n_valid > 1, (n_valid - 1) / n_valid * sumsq, np.nan)
    return mean, var, n_valid


mean_u, var_u, n_valid_u = jackknife_mean_var(tau_u_mat)
mean_w, var_w, n_valid_w = jackknife_mean_var(tau_w_mat)

summary = pd.DataFrame({
    "gene": genes,
    "n_reps_tau_unweighted": n_valid_u,
    "tau_unweighted_jk_mean": mean_u,
    "tau_unweighted_jk_var": var_u,
    "n_reps_tau_free_cluster_weighted": n_valid_w,
    "tau_free_cluster_weighted_jk_mean": mean_w,
    "tau_free_cluster_weighted_jk_var": var_w,
})
summary.to_csv(OUT_SUMMARY_CSV, index=False)

long_df = pd.DataFrame({
    "gene": np.tile(genes, n_tissues),
    "dropped_tissue": np.repeat(all_tissues, n_genes),
    "tau_unweighted_jk": tau_u_mat.T.flatten(),
    "tau_free_cluster_weighted_jk": tau_w_mat.T.flatten(),
})
long_df.to_csv(OUT_LONG_CSV, index=False)

pd.DataFrame(rho_log).to_csv(OUT_SUMMARY_CSV.replace(".csv", "_replicate_log.csv"), index=False)
print(f"wrote {OUT_SUMMARY_CSV}, {OUT_LONG_CSV}", file=sys.stderr)
