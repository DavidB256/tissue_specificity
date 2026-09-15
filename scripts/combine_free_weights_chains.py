"""Pick the best of N independent free_weights/free_cluster_weights DE chains (see
optimize_free_weights.py and the free_weights_chain snakemake rule) by final Spearman
rho, copy that chain's weights/tau (and, in cluster MODE, per-cluster weights) as the
canonical output, and write a summary of every chain's result for transparency."""
import sys
import shutil
import pandas as pd

MODE = sys.argv[1]
assert MODE in ("tissue", "cluster"), f"MODE must be tissue or cluster, got {MODE!r}"
CHAINS_DIR = sys.argv[2]
N_CHAINS = int(sys.argv[3])
OUT_WEIGHTS_CSV = sys.argv[4]
OUT_CLUSTER_WEIGHTS_CSV = sys.argv[5]  # "-" if MODE == "tissue"
OUT_TAU_CSV = sys.argv[6]
OUT_SUMMARY_CSV = sys.argv[7]

rows = []
for i in range(N_CHAINS):
    best = pd.read_csv(f"{CHAINS_DIR}/chain_{i}/best.csv").iloc[0]
    rows.append({"chain": i, "seed": int(best["seed"]), "best_rho": best["best_rho"],
                 "n_generations": int(best["n_generations"])})

summary = pd.DataFrame(rows).sort_values("best_rho", ascending=False)
summary.to_csv(OUT_SUMMARY_CSV, index=False)

winner = int(summary.iloc[0]["chain"])
shutil.copy(f"{CHAINS_DIR}/chain_{winner}/weights.csv", OUT_WEIGHTS_CSV)
shutil.copy(f"{CHAINS_DIR}/chain_{winner}/weighted_tau.csv", OUT_TAU_CSV)
if MODE == "cluster":
    shutil.copy(f"{CHAINS_DIR}/chain_{winner}/cluster_weights.csv", OUT_CLUSTER_WEIGHTS_CSV)

print(f"winner: chain={winner} rho={summary.iloc[0]['best_rho']:.6f} "
      f"(range across {N_CHAINS} chains: {summary['best_rho'].min():.6f} - {summary['best_rho'].max():.6f})",
      file=sys.stderr)
