""""free_weights": optimize one weight per tissue (68-dim, bounded [0.01, 1]) directly --
not derived from any GLS or clustering formula -- to maximize the Spearman rank
correlation between weighted tau on the full (68-tissue) GCT and the fixed reference
tau_unweighted on the onebrain (56-tissue, brain collapsed) GCT (already computed, see
compute_weighted_tau.py). L-BFGS is not viable here: tau_i(w) is a smooth ratio of two
linear functions of w (argmax per gene depends only on raw TPM, never on w), but Spearman
rho only depends on gene RANK order, which is piecewise-constant in w -- a small step in
w almost never swaps two genes' relative tau, so the gradient is exactly zero almost
everywhere, with jumps only at measure-zero rank-crossing boundaries. Uses
scipy.optimize.differential_evolution (derivative-free, global; see
optimize_cluster_weight_params.py for the same reasoning in the 2-D case) instead. One
run = one chain (SEED); run many chains in parallel (see the free_weights_chain
snakemake rule, wildcard-parallelized) and keep the best, since DE in 68 dimensions is
not guaranteed to find the same optimum from every random start."""
import sys
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from scipy.optimize import differential_evolution
from gls_tau import weighted_tau, floor_tpm

FULL_GCT = sys.argv[1]
ONEBRAIN_TAU_CSV = sys.argv[2]
OUT_WEIGHTS_CSV = sys.argv[3]
OUT_TAU_CSV = sys.argv[4]
OUT_LOG_CSV = sys.argv[5]
OUT_BEST_CSV = sys.argv[6]
SEED = int(sys.argv[7]) if len(sys.argv) > 7 else 0

MAXITER = 300
POPSIZE = 15

df = pd.read_csv(FULL_GCT, sep="\t", skiprows=2)
tissues = list(df.columns[2:])
n_tissues = len(tissues)
genes = df["Name"].values
mat = floor_tpm(df[tissues].values.astype(float))

ref = pd.read_csv(ONEBRAIN_TAU_CSV)[["gene", "tau_unweighted"]].rename(columns={"tau_unweighted": "tau_ref"})
ref = pd.DataFrame({"gene": genes}).merge(ref, on="gene", how="left")
tau_ref = ref["tau_ref"].values

# valid mask is fixed for all w: full-side NaN only comes from an all-zero
# expression row (never from w), so precompute once and slice for speed.
valid = np.isfinite(tau_ref) & (mat.max(axis=1) > 0)
mat_valid = mat[valid]
tau_ref_valid = tau_ref[valid]
print(f"tissues={n_tissues} genes={len(genes)} valid_genes={valid.sum()}", file=sys.stderr)

# weighted_tau recomputes max/argmax/the full (1-mat/max) term every call even though
# none of those depend on w -- precompute once (10x speedup, verified bit-identical to
# gls_tau.weighted_tau) since DE needs hundreds of thousands of evaluations.
max_i_valid = mat_valid.max(axis=1)
argmax_idx_valid = mat_valid.argmax(axis=1)
term_valid = 1 - mat_valid / max_i_valid[:, None]

rank_ref = rankdata(tau_ref_valid)
rank_ref_c = rank_ref - rank_ref.mean()
denom_ref = np.sqrt((rank_ref_c ** 2).sum())


def neg_spearman(w):
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


bounds = [(0.01, 1.0)] * n_tissues
result = differential_evolution(neg_spearman, bounds=bounds, seed=SEED, maxiter=MAXITER,
                                 popsize=POPSIZE, tol=1e-8, polish=False, disp=False,
                                 callback=log_generation)

pd.DataFrame(log_rows).to_csv(OUT_LOG_CSV, index=False)

w_best = result.x
best_rho = -result.fun
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
      f"w_min={w_best.min():.4f} w_max={w_best.max():.4f}", file=sys.stderr)
