""""free_cluster_weights": dimensionality-reduced variant of "free_weights"
(optimize_free_weights.py). Instead of one free DE parameter per tissue (68-dim for
GTEx), one free parameter is optimized per hierarchical-clustering cluster (see
CLUSTERS_CSV, e.g. cluster_weights/full/tissue_clusters.csv at DIST_THRESHOLD=0.1 --
15 clusters for GTEx's 68 tissues, 10 for Baron's 14 cell types) and broadcast to every
tissue in that cluster: w_tissue = w_cluster[cluster_id(tissue)]. Same objective as
optimize_free_weights.py -- maximize the Spearman rank correlation between
weighted_tau(full-resolution GCT, w) and the fixed reduced-resolution reference
tau_unweighted -- so the two methods are directly comparable: does clustering the
weight-choosing problem down to n_clusters dimensions cost much rho relative to the
full n_tissues-dimensional free_weights fit? Output weights.csv is broadcast to
per-tissue (same "tissue","weight" schema as optimize_free_weights.py) so every
downstream free_weights plotting/comparison script works unmodified; the raw
per-cluster fit is written separately to OUT_CLUSTER_WEIGHTS_CSV."""
import sys
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from scipy.optimize import differential_evolution
from gls_tau import weighted_tau, floor_tpm

FULL_GCT = sys.argv[1]
CLUSTERS_CSV = sys.argv[2]
REF_TAU_CSV = sys.argv[3]
OUT_WEIGHTS_CSV = sys.argv[4]
OUT_CLUSTER_WEIGHTS_CSV = sys.argv[5]
OUT_TAU_CSV = sys.argv[6]
OUT_LOG_CSV = sys.argv[7]
OUT_BEST_CSV = sys.argv[8]
SEED = int(sys.argv[9]) if len(sys.argv) > 9 else 0

MAXITER = 300
POPSIZE = 15

df = pd.read_csv(FULL_GCT, sep="\t", skiprows=2)
tissues = list(df.columns[2:])
n_tissues = len(tissues)
genes = df["Name"].values
mat = floor_tpm(df[tissues].values.astype(float))

clusters = pd.read_csv(CLUSTERS_CSV).set_index("tissue").loc[tissues]
cluster_ids = sorted(clusters["cluster_id"].unique())
cluster_of_tissue = clusters["cluster_id"].map({c: i for i, c in enumerate(cluster_ids)}).values
n_clusters = len(cluster_ids)

ref = pd.read_csv(REF_TAU_CSV)[["gene", "tau_unweighted"]].rename(columns={"tau_unweighted": "tau_ref"})
ref = pd.DataFrame({"gene": genes}).merge(ref, on="gene", how="left")
tau_ref = ref["tau_ref"].values

# valid mask is fixed for all w: full-side NaN only comes from an all-zero
# expression row (never from w), so precompute once and slice for speed.
valid = np.isfinite(tau_ref) & (mat.max(axis=1) > 0)
mat_valid = mat[valid]
tau_ref_valid = tau_ref[valid]
print(f"tissues={n_tissues} clusters={n_clusters} genes={len(genes)} valid_genes={valid.sum()}", file=sys.stderr)

# weighted_tau recomputes max/argmax/the full (1-mat/max) term every call even though
# none of those depend on w -- precompute once (10x speedup, verified bit-identical to
# gls_tau.weighted_tau) since DE needs hundreds of thousands of evaluations.
max_i_valid = mat_valid.max(axis=1)
argmax_idx_valid = mat_valid.argmax(axis=1)
term_valid = 1 - mat_valid / max_i_valid[:, None]

rank_ref = rankdata(tau_ref_valid)
rank_ref_c = rank_ref - rank_ref.mean()
denom_ref = np.sqrt((rank_ref_c ** 2).sum())


def neg_spearman(w_cluster):
    w = w_cluster[cluster_of_tissue]  # broadcast cluster weight -> every tissue in it
    tau = term_valid @ w / (w.sum() - w[argmax_idx_valid])
    rank_tau = rankdata(tau)
    rank_tau_c = rank_tau - rank_tau.mean()
    rho = (rank_tau_c @ rank_ref_c) / (np.sqrt((rank_tau_c ** 2).sum()) * denom_ref)
    return -rho


log_rows = []


def log_generation(xk, convergence):
    rho = -neg_spearman(xk)
    log_rows.append({"generation": len(log_rows), "best_rho": rho, "convergence": convergence})
    print(f"gen {len(log_rows)}: best_rho={rho:.6f} convergence={convergence:.4g}", file=sys.stderr)


bounds = [(0.01, 1.0)] * n_clusters
result = differential_evolution(neg_spearman, bounds=bounds, seed=SEED, maxiter=MAXITER,
                                 popsize=POPSIZE, tol=1e-8, polish=False, disp=False,
                                 callback=log_generation)

pd.DataFrame(log_rows).to_csv(OUT_LOG_CSV, index=False)

w_cluster_best = result.x
best_rho = -result.fun
pd.DataFrame({"cluster_id": cluster_ids, "weight": w_cluster_best}).to_csv(OUT_CLUSTER_WEIGHTS_CSV, index=False)

w_best = w_cluster_best[cluster_of_tissue]
pd.DataFrame({"tissue": tissues, "weight": w_best}).to_csv(OUT_WEIGHTS_CSV, index=False)

tau_free = weighted_tau(mat, w_best)
out = pd.DataFrame({
    "gene": genes,
    "tau_unweighted": weighted_tau(mat, np.ones(n_tissues)),
    "tau_free_weighted": tau_free,
})
out.to_csv(OUT_TAU_CSV, index=False)

pd.DataFrame([{"seed": SEED, "best_rho": best_rho, "n_generations": len(log_rows)}]).to_csv(OUT_BEST_CSV, index=False)

print(f"seed={SEED} best: rho={best_rho:.6f} n_generations={len(log_rows)} "
      f"n_clusters={n_clusters} w_min={w_best.min():.4f} w_max={w_best.max():.4f}", file=sys.stderr)
