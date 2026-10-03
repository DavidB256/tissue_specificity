"""2-subfig scatter: (i) gene median_log_tpm (log10(TPM+1), median across
tissues, no flooring) vs. gene tau_unweighted, (ii) same for transcripts
(median_log_tpm computed from the protein-coding GCT) vs. transcript
tau_unweighted. Points colored by n_transcripts (log scale) -- the gene's own
transcript count in panel (i), the gene-of-that-transcript's transcript count
in panel (ii). Genes/transcripts with TPM<1 in every tissue are excluded --
tau_unweighted is NaN for exactly that all-zero-row set (gtex_v11_tx floors
TPM<1 to 0 upstream), so the existing NaN-tau dropna already enforces it; this
was verified by direct max-TPM-across-tissues count match. Red OLS fit
y = b0 + b1*x (statsmodels) on each panel, equation + 95% CIs for b0 and b1 in
the legend."""
import sys
import numpy as np
import pandas as pd
import statsmodels.api as sm
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

def ols_fit_line(ax, x, y):
    model = sm.OLS(y, sm.add_constant(x)).fit()
    b0, b1 = model.params
    (b0_lo, b0_hi), (b1_lo, b1_hi) = model.conf_int(alpha=0.05)
    x_line = np.array([x.min(), x.max()])
    label = (
        f"y = {b0:.3f} + {b1:.3f} x\n"
        f"95% CI $\\beta_0$: [{b0_lo:.3f}, {b0_hi:.3f}]\n"
        f"95% CI $\\beta_1$: [{b1_lo:.3f}, {b1_hi:.3f}]"
    )
    ax.plot(x_line, b0 + b1 * x_line, color="red", linewidth=2, label=label)
    ax.legend(loc="lower right", fontsize=9)
    return b0, b1, b0_lo, b0_hi, b1_lo, b1_hi


fig, axes = plt.subplots(1, 2, figsize=(14, 7))

sc0 = axes[0].scatter(gene_df["median_log_tpm"], gene_df["tau_unweighted"], c=gene_df["n_transcripts"],
                       cmap="viridis", s=10, alpha=0.3, norm=LogNorm())
stats0 = ols_fit_line(axes[0], gene_df["median_log_tpm"].to_numpy(), gene_df["tau_unweighted"].to_numpy())
axes[0].set_xlabel("gene median_log_tpm (log10 TPM+1)")
axes[0].set_ylabel("gene tau_unweighted")
axes[0].set_title(f"gene log TPM vs. gene tau (n={len(gene_df)} genes)")
fig.colorbar(sc0, ax=axes[0], label="n_transcripts (log scale)")

sc1 = axes[1].scatter(tx_df["median_log_tpm"], tx_df["tau_tx"], c=tx_df["n_transcripts"],
                       cmap="viridis", s=10, alpha=0.1, norm=LogNorm())
stats1 = ols_fit_line(axes[1], tx_df["median_log_tpm"].to_numpy(), tx_df["tau_tx"].to_numpy())
axes[1].set_xlabel("transcript median_log_tpm (log10 TPM+1)")
axes[1].set_ylabel("transcript tau_unweighted")
axes[1].set_title(f"tx log TPM vs. tx tau (n={len(tx_df)} transcripts)")
fig.colorbar(sc1, ax=axes[1], label="n_transcripts of tx's gene (log scale)")

fig.suptitle(f"{TXS_SOURCE}: log10(TPM+1) vs. tau (TPM<1-in-all-tissues excluded)")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)

for label, n, (b0, b1, b0_lo, b0_hi, b1_lo, b1_hi) in [
    ("gene", len(gene_df), stats0), ("tx", len(tx_df), stats1)
]:
    print(f"{label}: n={n} b0={b0:.4f} [{b0_lo:.4f},{b0_hi:.4f}] b1={b1:.4f} [{b1_lo:.4f},{b1_hi:.4f}]",
          file=sys.stderr)
print(f"-> {OUT_PNG}", file=sys.stderr)
