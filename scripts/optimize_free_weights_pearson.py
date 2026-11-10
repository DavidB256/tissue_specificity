"""free_cluster_weights, Pearson variant: maximize Pearson r (not Spearman rho) between
weighted tau on the full-resolution GCT and a fixed reference tau (see
optimize_free_weights.py for the problem, cluster parametrization and arguments; same
CLI and output files, plus runtime/eval counts in OUT_BEST_CSV). r is smooth in w, so
the gradient comes from torch autograd (float64, CPU) and the optimizer is scipy
L-BFGS-B over bounded cluster weights. One run = one chain = N_STARTS L-BFGS starts,
keeping the best; the first start of chain SEED==0 is the naive guess
1/cluster_size, every other start is log-uniform in the bounds."""
import sys
import time
import numpy as np
import pandas as pd
import torch
from scipy.stats import rankdata, spearmanr
from scipy.optimize import minimize
from scipy.cluster.hierarchy import fcluster
from gls_tau import weighted_tau, floor_tpm, log_tpm

FULL_GCT = sys.argv[1]
REF_TAU_CSV = sys.argv[2]
CLUSTERS_CSV = sys.argv[3]  # "-" -> cluster LINKAGE_CSV at DIST_THRESHOLD
OUT_WEIGHTS_CSV = sys.argv[4]
OUT_CLUSTER_WEIGHTS_CSV = sys.argv[5]
OUT_TAU_CSV = sys.argv[6]
OUT_LOG_CSV = sys.argv[7]
OUT_BEST_CSV = sys.argv[8]
SEED = int(sys.argv[9]) if len(sys.argv) > 9 else 0
DIST_THRESHOLD = float(sys.argv[10]) if len(sys.argv) > 10 else None
LINKAGE_CSV = sys.argv[11] if len(sys.argv) > 11 else None

N_STARTS = 20
LB, UB = 0.01, 1.0
t0 = time.perf_counter()

df = pd.read_csv(FULL_GCT, sep="\t", skiprows=2)
tissues = list(df.columns[2:])
n_tissues = len(tissues)
genes = df["Name"].values
mat = floor_tpm(df[tissues].values.astype(float))

if CLUSTERS_CSV == "-":
    Z = pd.read_csv(LINKAGE_CSV)[["cluster1", "cluster2", "distance", "n_obs"]].values
    raw_ids = fcluster(Z, t=DIST_THRESHOLD, criterion="distance")
else:
    raw_ids = pd.read_csv(CLUSTERS_CSV).set_index("tissue").loc[tissues]["cluster_id"].values
cluster_ids = sorted(set(raw_ids))
cluster_of_tissue = np.array([cluster_ids.index(c) for c in raw_ids])
n_dims = len(cluster_ids)

ref = pd.read_csv(REF_TAU_CSV)[["gene", "tau_unweighted"]].rename(columns={"tau_unweighted": "tau_ref"})
tau_ref = pd.DataFrame({"gene": genes}).merge(ref, on="gene", how="left")["tau_ref"].values
valid = np.isfinite(tau_ref) & (mat.max(axis=1) > 0)
mat_valid = mat[valid]
tau_ref_valid = tau_ref[valid]
print(f"tissues={n_tissues} dims={n_dims} genes={len(genes)} valid_genes={valid.sum()}", file=sys.stderr)

# tau_g(w) = term_g . w / (sum(w) - w[argmax_g]); term/argmax never depend on w
log_valid = log_tpm(mat_valid)
argmax_idx = torch.from_numpy(log_valid.argmax(axis=1))
term = torch.from_numpy(1 - log_valid / log_valid.max(axis=1)[:, None])
cluster_idx = torch.from_numpy(cluster_of_tissue)
rc = torch.from_numpy(tau_ref_valid - tau_ref_valid.mean())
rc_norm = rc.norm()


def neg_r_and_grad(x):
    xt = torch.tensor(x, dtype=torch.float64, requires_grad=True)
    w = xt[cluster_idx]
    tau = term @ w / (w.sum() - w[argmax_idx])
    tc = tau - tau.mean()
    r = (tc @ rc) / (tc.norm() * rc_norm)
    (-r).backward()
    return -r.item(), xt.grad.numpy().copy()


# gradient gate: central finite difference vs autograd
x_chk = np.random.default_rng(SEED).uniform(LB, UB, n_dims)
g_ad = neg_r_and_grad(x_chk)[1]
g_fd = np.array([(neg_r_and_grad(x_chk + h)[0] - neg_r_and_grad(x_chk - h)[0]) / (2e-5)
                 for h in 1e-5 * np.eye(n_dims)])
rel = np.linalg.norm(g_ad - g_fd) / max(np.linalg.norm(g_fd), 1e-12)
print(f"grad check: rel err autograd vs finite diff = {rel:.2e}", file=sys.stderr)
assert rel < 1e-4

rng = np.random.default_rng(SEED)
starts = [np.exp(rng.uniform(np.log(LB), np.log(UB), n_dims)) for _ in range(N_STARTS)]
if SEED == 0:
    starts[0] = np.clip(1.0 / np.bincount(cluster_of_tissue, minlength=n_dims), LB, UB)

log_rows, best = [], None
for k, x0 in enumerate(starts):
    res = minimize(neg_r_and_grad, x0, jac=True, method="L-BFGS-B", bounds=[(LB, UB)] * n_dims,
                   options=dict(maxiter=1000, ftol=1e-13, gtol=1e-10))
    log_rows.append({"start": k, "r": -res.fun, "nfev": res.nfev, "nit": res.nit})
    print(f"start {k}: r={-res.fun:.8f} nfev={res.nfev}", file=sys.stderr)
    if best is None or res.fun < best.fun:
        best = res
pd.DataFrame(log_rows).to_csv(OUT_LOG_CSV, index=False)

x_best = best.x
w_best = x_best[cluster_of_tissue]
pd.DataFrame({"cluster_id": cluster_ids, "weight": x_best}).to_csv(OUT_CLUSTER_WEIGHTS_CSV, index=False)
pd.DataFrame({"tissue": tissues, "weight": w_best}).to_csv(OUT_WEIGHTS_CSV, index=False)

tau_free = weighted_tau(mat, w_best)
pd.DataFrame({
    "gene": genes,
    "tau_unweighted": weighted_tau(mat, np.ones(n_tissues)),
    "tau_free_weighted": tau_free,
}).to_csv(OUT_TAU_CSV, index=False)

rho = spearmanr(tau_free[valid], tau_ref_valid)[0]
runtime = time.perf_counter() - t0
pd.DataFrame([{"seed": SEED, "best_r": -best.fun, "rho_at_best": rho, "n_starts": N_STARTS,
               "n_evals": sum(r["nfev"] for r in log_rows), "runtime_s": runtime}]).to_csv(OUT_BEST_CSV, index=False)
print(f"seed={SEED} best: r={-best.fun:.8f} rho={rho:.6f} runtime={runtime:.1f}s "
      f"w_min={w_best.min():.4f} w_max={w_best.max():.4f}", file=sys.stderr)
