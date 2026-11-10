"""Pick the best of N independent free_cluster_weights DE chains (see
optimize_free_weights.py and the optimize_free_cluster_weights_chain snakemake rule) by final Spearman
rho, copy that chain's weights/tau/per-cluster weights as the
canonical output, and write a summary of every chain's result for transparency."""
import sys
import shutil
import pandas as pd

CHAINS_DIR = sys.argv[1]
N_CHAINS = int(sys.argv[2])
OUT_WEIGHTS_CSV = sys.argv[3]
OUT_CLUSTER_WEIGHTS_CSV = sys.argv[4]
OUT_TAU_CSV = sys.argv[5]
OUT_SUMMARY_CSV = sys.argv[6]
METRIC = sys.argv[7] if len(sys.argv) > 7 else "best_rho"  # best.csv column to maximize

rows = []
for i in range(N_CHAINS):
    best = pd.read_csv(f"{CHAINS_DIR}/chain_{i}/best.csv").iloc[0]
    rows.append({"chain": i, **best.to_dict()})

summary = pd.DataFrame(rows).sort_values(METRIC, ascending=False)
summary.to_csv(OUT_SUMMARY_CSV, index=False)

winner = int(summary.iloc[0]["chain"])
shutil.copy(f"{CHAINS_DIR}/chain_{winner}/weights.csv", OUT_WEIGHTS_CSV)
shutil.copy(f"{CHAINS_DIR}/chain_{winner}/weighted_tau.csv", OUT_TAU_CSV)
shutil.copy(f"{CHAINS_DIR}/chain_{winner}/cluster_weights.csv", OUT_CLUSTER_WEIGHTS_CSV)

print(f"winner: chain={winner} {METRIC}={summary.iloc[0][METRIC]:.6f} "
      f"(range across {N_CHAINS} chains: {summary[METRIC].min():.6f} - {summary[METRIC].max():.6f})",
      file=sys.stderr)
