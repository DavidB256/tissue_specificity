"""2-subfig scatter: (i) tx-gene winner distance vs. gene tau_unweighted, (ii)
tx-gene winner distance vs. transcript tau_unweighted. Distance is
1 - Pearson r (log10(TPM+1), TPM<1->0) between a transcript's winning tissue
and its gene's winning tissue (compute_tx_gene_winner_distance.py)."""
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IN_DIST_CSV = sys.argv[1]
IN_GENE_TAU_CSV = sys.argv[2]
IN_TX_TAU_CSV = sys.argv[3]
OUT_PNG = sys.argv[4]
TXS_SOURCE = sys.argv[5]

dist = pd.read_csv(IN_DIST_CSV)[["transcript_id", "gene_id", "dist"]]
gene_tau = pd.read_csv(IN_GENE_TAU_CSV)[["gene_id", "tau_unweighted"]].rename(columns={"tau_unweighted": "tau_gene"})
tx_tau = pd.read_csv(IN_TX_TAU_CSV).rename(columns={"gene": "transcript_id", "tau_unweighted": "tau_tx"})

df = dist.merge(gene_tau, on="gene_id", how="left").merge(tx_tau, on="transcript_id", how="left")

fig, axes = plt.subplots(1, 2, figsize=(14, 7))

d0 = df.dropna(subset=["tau_gene", "dist"])
axes[0].scatter(d0["tau_gene"], d0["dist"], s=8, alpha=0.1, color="steelblue")
axes[0].set_xlabel("gene tau_unweighted")
axes[0].set_ylabel("tx-gene winner distance (1 - Pearson r)")
axes[0].set_title(f"distance vs. gene tau (n={len(d0)})")

d1 = df.dropna(subset=["tau_tx", "dist"])
axes[1].scatter(d1["tau_tx"], d1["dist"], s=8, alpha=0.1, color="steelblue")
axes[1].set_xlabel("transcript tau_unweighted")
axes[1].set_ylabel("tx-gene winner distance (1 - Pearson r)")
axes[1].set_title(f"distance vs. transcript tau (n={len(d1)})")

fig.suptitle(f"{TXS_SOURCE}: tx-gene winning-tissue distance vs. tau")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)

print(f"n_gene_tau={len(d0)} n_tx_tau={len(d1)} -> {OUT_PNG}", file=sys.stderr)
