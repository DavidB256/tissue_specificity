"""10x10 grid of KDE plots, one per gene, for the 100 protein-coding gtex_txs genes
with the most transcripts (multitranscript_variance.csv's n_transcripts, descending).
Each panel shows the distribution of tau_unweighted (compute_tau_unweighted.py) among
that gene's transcripts, labeled with the gene's name (CHESS gene_name, not gene_id)
and n_txs. Genes with fewer than 2 transcripts having a defined (non-NaN) tau fall
back to a rug of the raw values (too little data for a KDE)."""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

N_GENES = 100
GRID = 10

IN_TAU_CSV = sys.argv[1]
IN_META_CSV = sys.argv[2]
IN_MULTITX_CSV = sys.argv[3]
OUT_PNG = sys.argv[4]

tau = pd.read_csv(IN_TAU_CSV).rename(columns={"gene": "transcript_id"})
meta = pd.read_csv(IN_META_CSV)[["transcript_id", "gene_id", "gene_name"]]
tau = tau.merge(meta, on="transcript_id", how="left")

multitx = pd.read_csv(IN_MULTITX_CSV)
top_genes = multitx.sort_values("n_transcripts", ascending=False).head(N_GENES)

gene_name = meta.drop_duplicates("gene_id").set_index("gene_id")["gene_name"]

fig, axes = plt.subplots(GRID, GRID, figsize=(28, 28))
for ax, row in zip(axes.flat, top_genes.itertuples()):
    vals = tau.loc[tau["gene_id"] == row.gene_id, "tau_unweighted"].dropna().to_numpy()
    name = gene_name.get(row.gene_id, row.gene_id)
    if len(vals) >= 2 and vals.std() > 0:
        sns.kdeplot(vals, ax=ax, fill=True, color="steelblue")
    else:
        ax.plot(vals, np.zeros_like(vals), "|", color="steelblue", markersize=20, markeredgewidth=2)
    ax.set_title(f"{name} (n_txs={row.n_transcripts})", fontsize=9)
    ax.set_xlim(0, 1)
    ax.set_xlabel("")
    ax.set_ylabel("")

fig.suptitle(f"gtex_txs: tau distribution among transcripts, top {N_GENES} genes by n_transcripts", y=1.0)
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)

print(f"genes={len(top_genes)} -> {OUT_PNG}", file=sys.stderr)
