"""2-subfig scatter: (i) gene median_log_tpm (log10(TPM+1), median across
tissues, no flooring) vs. gene tau_unweighted, (ii) same for transcripts
(median_log_tpm computed from the protein-coding GCT) vs. transcript
tau_unweighted. Points colored by n_transcripts (log scale) -- the gene's own
transcript count in panel (i), the gene-of-that-transcript's transcript count
in panel (ii)."""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

IN_GENE_TAU_CSV = sys.argv[1]
IN_TX_TAU_CSV = sys.argv[2]
IN_GCT = sys.argv[3]
OUT_PNG = sys.argv[4]
TXS_SOURCE = sys.argv[5]

gct = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(gct.columns[2:])
tx_mat = gct[tissues].to_numpy()
tx_log_tpm = pd.DataFrame({
    "transcript_id": gct["Name"],
    "gene_id": gct["Description"],
    "median_log_tpm": np.median(np.log10(tx_mat + 1.0), axis=1),
})
n_transcripts = gct.groupby("Description").size().rename("n_transcripts")

gene_tau = pd.read_csv(IN_GENE_TAU_CSV)[["gene_id", "tau_unweighted", "median_log_tpm"]]
gene_df = gene_tau.merge(n_transcripts, left_on="gene_id", right_index=True, how="left").dropna(subset=["tau_unweighted"])

tx_tau = pd.read_csv(IN_TX_TAU_CSV).rename(columns={"gene": "transcript_id", "tau_unweighted": "tau_tx"})
tx_df = tx_log_tpm.merge(tx_tau, on="transcript_id", how="left")
tx_df = tx_df.merge(n_transcripts, left_on="gene_id", right_index=True, how="left").dropna(subset=["tau_tx"])

fig, axes = plt.subplots(1, 2, figsize=(14, 7))

sc0 = axes[0].scatter(gene_df["median_log_tpm"], gene_df["tau_unweighted"], c=gene_df["n_transcripts"],
                       cmap="viridis", s=10, alpha=0.3, norm=LogNorm())
axes[0].set_xlabel("gene median_log_tpm (log10 TPM+1)")
axes[0].set_ylabel("gene tau_unweighted")
axes[0].set_title(f"gene log TPM vs. gene tau (n={len(gene_df)} genes)")
fig.colorbar(sc0, ax=axes[0], label="n_transcripts (log scale)")

sc1 = axes[1].scatter(tx_df["median_log_tpm"], tx_df["tau_tx"], c=tx_df["n_transcripts"],
                       cmap="viridis", s=10, alpha=0.1, norm=LogNorm())
axes[1].set_xlabel("transcript median_log_tpm (log10 TPM+1)")
axes[1].set_ylabel("transcript tau_unweighted")
axes[1].set_title(f"tx log TPM vs. tx tau (n={len(tx_df)} transcripts)")
fig.colorbar(sc1, ax=axes[1], label="n_transcripts of tx's gene (log scale)")

fig.suptitle(f"{TXS_SOURCE}: log10(TPM+1) vs. tau")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)

print(f"genes={len(gene_df)} transcripts={len(tx_df)} -> {OUT_PNG}", file=sys.stderr)
