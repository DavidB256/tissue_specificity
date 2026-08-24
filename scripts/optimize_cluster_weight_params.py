"""Optimize (DIST_THRESHOLD, WEIGHT_CAP) for the diameter-scaled cluster weight
weight_diameter = (1-d)/cluster_size + WEIGHT_CAP*d (see compute_tau_cluster_weighted.py),
maximizing the Spearman rank correlation of tau_weighted_cluster_diameter between the
default (all 68 tissues) and onebrain (56 tissues, brain collapsed) tissue sets --
i.e. the params that make cluster-weighted tau most rank-consistent across tissue
granularities (research goal: tau shouldn't depend on how finely tissues are split).
Uses scipy.optimize.differential_evolution (derivative-free, global) since the objective
is piecewise-constant in DIST_THRESHOLD (clustering only changes at merge-height
crossings) and only piecewise-linear in WEIGHT_CAP given fixed clusters -- not smooth
enough for gradient methods. Every evaluated (threshold, cap, rho) triple is logged to
OUT_TRAJECTORY_CSV; the best is written to OUT_BEST_CSV."""
import sys
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster
from scipy.stats import spearmanr
from scipy.optimize import differential_evolution
from gls_tau import weighted_tau, floor_tpm

DEFAULT_GCT = sys.argv[1]
DEFAULT_LINKAGE = sys.argv[2]
DEFAULT_CORR = sys.argv[3]
ONEBRAIN_GCT = sys.argv[4]
ONEBRAIN_LINKAGE = sys.argv[5]
ONEBRAIN_CORR = sys.argv[6]
OUT_TRAJECTORY_CSV = sys.argv[7]
OUT_BEST_CSV = sys.argv[8]

BOUNDS = [(0.01, 0.5), (0.1, 0.9)]  # (dist_threshold, weight_cap)
SEED = 0


def load_leg(gct_path, linkage_path, corr_path):
    df = pd.read_csv(gct_path, sep="\t", skiprows=2)
    tissues = list(df.columns[2:])
    genes = df["Name"].values
    mat = floor_tpm(df[tissues].values.astype(float))
    Z = pd.read_csv(linkage_path)[["cluster1", "cluster2", "distance", "n_obs"]].values
    R = pd.read_csv(corr_path, index_col=0).loc[tissues, tissues].values
    dist_mat = 1 - R
    return genes, mat, Z, dist_mat


genes_d, mat_d, Z_d, dist_d = load_leg(DEFAULT_GCT, DEFAULT_LINKAGE, DEFAULT_CORR)
genes_o, mat_o, Z_o, dist_o = load_leg(ONEBRAIN_GCT, ONEBRAIN_LINKAGE, ONEBRAIN_CORR)
assert (genes_d == genes_o).all(), "gene order mismatch between default and onebrain GCTs"


def cluster_weights(Z, dist_mat, n_tissues, threshold, cap):
    cluster_id = fcluster(Z, t=threshold, criterion="distance")
    cluster_size = pd.Series(cluster_id).value_counts()
    size_arr = np.array([cluster_size[c] for c in cluster_id])
    diameter = np.zeros(n_tissues)
    for c in cluster_size.index:
        idx = np.where(cluster_id == c)[0]
        if len(idx) > 1:
            diameter[idx] = dist_mat[np.ix_(idx, idx)].max()
    d = np.minimum(diameter / threshold, 1.0)
    return (1 - d) / size_arr + cap * d


trajectory = []


def neg_spearman(x):
    threshold, cap = x
    w_d = cluster_weights(Z_d, dist_d, mat_d.shape[1], threshold, cap)
    w_o = cluster_weights(Z_o, dist_o, mat_o.shape[1], threshold, cap)
    tau_d = weighted_tau(mat_d, w_d)
    tau_o = weighted_tau(mat_o, w_o)
    valid = np.isfinite(tau_d) & np.isfinite(tau_o)
    rho, _ = spearmanr(tau_d[valid], tau_o[valid])
    trajectory.append({"dist_threshold": threshold, "weight_cap": cap, "spearman_rho": rho, "n_genes": int(valid.sum())})
    return -rho


result = differential_evolution(neg_spearman, bounds=BOUNDS, seed=SEED, maxiter=60,
                                 popsize=15, tol=1e-6, polish=True, disp=False)

traj_df = pd.DataFrame(trajectory)
traj_df.to_csv(OUT_TRAJECTORY_CSV, index=False)

best_threshold, best_cap = result.x
best_rho = -result.fun
pd.DataFrame([{"dist_threshold": best_threshold, "weight_cap": best_cap, "spearman_rho": best_rho,
               "n_evaluations": len(trajectory)}]).to_csv(OUT_BEST_CSV, index=False)

print(f"best: threshold={best_threshold:.4f} cap={best_cap:.4f} rho={best_rho:.6f} "
      f"n_evaluations={len(trajectory)}", file=sys.stderr)
