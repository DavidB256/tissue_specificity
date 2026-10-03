"""3-subfig scatter: one point per gene, x=gene tau_unweighted
(compute_gene_tau_gtex_txs.py). Panel 1 (all genes) and panel 2
(multi-transcript genes only, >=2 protein-coding transcripts per
multitranscript_variance.csv): y=plain mean tau_unweighted across the gene's
transcripts (compute_tau_unweighted.py, NaN tx tau dropped before averaging).
Panel 3 (multi-transcript genes only): y=abundance-weighted mean tx tau --
each transcript's tau weighted by its mean TPM across tissues (from the
protein-coding GCT), e.g. a gene with tx1 at 3x tx2's expression gets
y=(3*tau1+tau2)/4. OLS fit y = b0 + b1*x (statsmodels) on each, equation +
95% CIs for b0 and b1 reported in the legend."""
import sys
import numpy as np
import pandas as pd
import statsmodels.api as sm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IN_TX_TAU_CSV = sys.argv[1]
IN_META_CSV = sys.argv[2]
IN_GENE_TAU_CSV = sys.argv[3]
IN_MULTITX_CSV = sys.argv[4]
IN_GCT = sys.argv[5]
OUT_PNG = sys.argv[6]
TXS_SOURCE = sys.argv[7]

tx_tau = pd.read_csv(IN_TX_TAU_CSV).rename(columns={"gene": "transcript_id", "tau_unweighted": "tau_tx"})
meta = pd.read_csv(IN_META_CSV)[["transcript_id", "gene_id"]]
gene_tau = pd.read_csv(IN_GENE_TAU_CSV)[["gene_id", "tau_unweighted"]].rename(columns={"tau_unweighted": "tau_gene"})
multitx_genes = pd.read_csv(IN_MULTITX_CSV)["gene_id"]

gct = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(gct.columns[2:])
abundance = pd.DataFrame({
    "transcript_id": gct["Name"],
    "abundance": gct[tissues].to_numpy().mean(axis=1),
})

tx = tx_tau.merge(meta, on="transcript_id", how="left").dropna(subset=["tau_tx"])
mean_tx_tau = tx.groupby("gene_id")["tau_tx"].mean().rename("mean_tau_tx").reset_index()

tx_w = tx.merge(abundance, on="transcript_id", how="left")
tx_w["contrib"] = tx_w["abundance"] * tx_w["tau_tx"]
weight_sum = tx_w.groupby("gene_id")["abundance"].sum()
contrib_sum = tx_w.groupby("gene_id")["contrib"].sum()
weighted_tx_tau = (contrib_sum / weight_sum.where(weight_sum > 0)).rename("weighted_tau_tx").reset_index()

df_all = gene_tau.merge(mean_tx_tau, on="gene_id", how="inner").dropna(subset=["tau_gene", "mean_tau_tx"])
df_multi = df_all[df_all["gene_id"].isin(multitx_genes)]
df_multi_weighted = (
    gene_tau.merge(weighted_tx_tau, on="gene_id", how="inner")
    .dropna(subset=["tau_gene", "weighted_tau_tx"])
)
df_multi_weighted = df_multi_weighted[df_multi_weighted["gene_id"].isin(multitx_genes)]


def fit_and_plot(ax, x, y, subtitle, ylabel):
    model = sm.OLS(y, sm.add_constant(x)).fit()
    b0, b1 = model.params
    (b0_lo, b0_hi), (b1_lo, b1_hi) = model.conf_int(alpha=0.05)

    ax.scatter(x, y, s=8, alpha=0.1, color="steelblue")
    x_line = np.array([x.min(), x.max()])
    label = (
        f"y = {b0:.3f} + {b1:.3f} x\n"
        f"95% CI $\\beta_0$: [{b0_lo:.3f}, {b0_hi:.3f}]\n"
        f"95% CI $\\beta_1$: [{b1_lo:.3f}, {b1_hi:.3f}]"
    )
    ax.plot(x_line, b0 + b1 * x_line, color="red", linewidth=2, label=label)
    ax.set_xlabel("gene tau_unweighted")
    ax.set_ylabel(ylabel)
    ax.set_title(f"{subtitle} (n={len(x)} genes)")
    ax.legend(loc="lower right", fontsize=9)
    return len(x), b0, b1, b0_lo, b0_hi, b1_lo, b1_hi


fig, axes = plt.subplots(1, 3, figsize=(21, 7))
stats_all = fit_and_plot(axes[0], df_all["tau_gene"].to_numpy(), df_all["mean_tau_tx"].to_numpy(),
                          "all genes, plain mean", "mean transcript tau_unweighted")
stats_multi = fit_and_plot(axes[1], df_multi["tau_gene"].to_numpy(), df_multi["mean_tau_tx"].to_numpy(),
                            "multi-transcript genes, plain mean", "mean transcript tau_unweighted")
stats_multi_w = fit_and_plot(axes[2], df_multi_weighted["tau_gene"].to_numpy(),
                              df_multi_weighted["weighted_tau_tx"].to_numpy(),
                              "multi-transcript genes, abundance-weighted",
                              "abundance-weighted mean transcript tau_unweighted")

fig.suptitle(f"{TXS_SOURCE}: gene tau vs. mean transcript tau")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)

for label, (n, b0, b1, b0_lo, b0_hi, b1_lo, b1_hi) in [
    ("all_plain", stats_all), ("multi_plain", stats_multi), ("multi_weighted", stats_multi_w)
]:
    print(f"{label}: n={n} b0={b0:.4f} [{b0_lo:.4f},{b0_hi:.4f}] b1={b1:.4f} [{b1_lo:.4f},{b1_hi:.4f}]",
          file=sys.stderr)
print(f"-> {OUT_PNG}", file=sys.stderr)
