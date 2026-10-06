"""Violin: x = whether EVERY transcript of a multi-tx gene (n_tx>=2) has
prob_tx_tau_gt_gene_weighted >= 0.95 (bootstrap, gtex_v11_tx_protein_coding),
y = gene_tau_weighted. Restricted to multi-tx genes (n_tx>=2) throughout --
single-tx genes trivially have tx tau == gene tau, so the grouping is
meaningless for them."""
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IN_BOOTSTRAP_TX_TAU_CSV = sys.argv[1]
IN_INFO_CONTENT_CSV = sys.argv[2]
OUT_PNG = sys.argv[3]

tx_tau = pd.read_csv(IN_BOOTSTRAP_TX_TAU_CSV)[["gene_id", "prob_tx_tau_gt_gene_weighted"]]
info = pd.read_csv(IN_INFO_CONTENT_CSV)[["gene_id", "n_tx", "gene_tau_weighted"]]

df = tx_tau.merge(info, on="gene_id", how="inner")
multi = df[df["n_tx"] >= 2]

all_sig = multi.groupby("gene_id")["prob_tx_tau_gt_gene_weighted"].apply(lambda g: (g >= 0.95).all())
all_sig.name = "all_tx_sig"

genes = info[info["n_tx"] >= 2][["gene_id", "gene_tau_weighted"]].merge(
    all_sig, left_on="gene_id", right_index=True, how="inner"
)

groups = [
    genes.loc[~genes["all_tx_sig"], "gene_tau_weighted"].dropna(),
    genes.loc[genes["all_tx_sig"], "gene_tau_weighted"].dropna(),
]
labels = [f"not all tx sig\n(n={len(groups[0])})", f"all tx sig\n(n={len(groups[1])})"]

fig, ax = plt.subplots(figsize=(7, 7))
parts = ax.violinplot(groups, positions=[1, 2], showmedians=True, widths=0.7)
for pc in parts["bodies"]:
    pc.set_facecolor("steelblue")
    pc.set_alpha(0.5)
ax.set_xticks([1, 2])
ax.set_xticklabels(labels)
ax.set_ylabel("gene tau_weighted")
ax.set_title("gtex_v11_tx_protein_coding: gene tau_weighted by all-tx-significance\n"
             "(every tx's prob_tx_tau_gt_gene_weighted >= 0.95), multi-tx genes only")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)

print(f"not_all_sig: n={len(groups[0])} median={groups[0].median():.4f} | "
      f"all_sig: n={len(groups[1])} median={groups[1].median():.4f} -> {OUT_PNG}", file=sys.stderr)
