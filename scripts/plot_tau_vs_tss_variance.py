"""Gene tau_weighted vs. variance in TSS location across a gene's transcripts.
TSS variance computed only over transcripts in transcript_meta (i.e. the
quantified gtex_v11_tx_protein_coding set), not every GENCODE transcript of
the gene. Genes with <2 such transcripts dropped (variance undefined)."""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

IN_TSS_CSV = sys.argv[1]
IN_TRANSCRIPT_META_CSV = sys.argv[2]
IN_INFO_CONTENT_CSV = sys.argv[3]
OUT_PNG = sys.argv[4]

tss = pd.read_csv(IN_TSS_CSV)[["transcript_id", "tss"]]
meta = pd.read_csv(IN_TRANSCRIPT_META_CSV)[["transcript_id", "gene_id"]]
info = pd.read_csv(IN_INFO_CONTENT_CSV)[["gene_id", "n_tx", "gene_tau_weighted"]]

df = tss.merge(meta, on="transcript_id", how="inner")
tss_var = df.groupby("gene_id")["tss"].agg(tss_variance="var", n_tss="count")
tss_var = tss_var[tss_var["n_tss"] >= 2]

gene = tss_var.join(info.set_index("gene_id"), how="inner").dropna(subset=["tss_variance", "gene_tau_weighted"])
gene["log10_tss_variance"] = np.log10(gene["tss_variance"] + 1)

rho, pval = spearmanr(gene["gene_tau_weighted"], gene["tss_variance"])

fig, ax = plt.subplots(figsize=(8, 7))
ax.scatter(gene["gene_tau_weighted"], gene["log10_tss_variance"], s=5, alpha=0.3)
ax.set_xlabel("gene tau_weighted")
ax.set_ylabel("log10(Var(TSS) + 1)  [bp^2]")
ax.set_title(
    f"gtex_v11_tx_protein_coding: gene tau_weighted vs. TSS variance across transcripts\n"
    f"n_genes={len(gene)} (n_tx>=2), spearman rho={rho:.3f} (p={pval:.2g})"
)
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)

print(f"n_genes={len(gene)} spearman_rho={rho:.4f} pval={pval:.4g} -> {OUT_PNG}", file=sys.stderr)
