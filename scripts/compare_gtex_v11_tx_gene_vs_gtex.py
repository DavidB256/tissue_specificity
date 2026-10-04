"""One-off check: are per-gene tau values from gtex_v11_tx_gene (RSEM v11
transcript TPMs summed to gene level) nearly identical to gtex_protein_coding
(RNA-SeQC's own gene-level TPM, the main pipeline's GTEx source)? Two
independent GTEx quantifications of the same underlying samples -- if the
pipeline is sound, gene-level tau should agree closely despite different
quantifiers (RSEM vs RNA-SeQC) and sample sets (v11 vs the main GCT's
release). Base (version-stripped) gene_id match. Not wired into snakefile --
ad hoc verification, not a pipeline output."""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IN_V11_TX_GENE_TAU = sys.argv[1]
IN_GTEX_TAU = sys.argv[2]
OUT_PNG = sys.argv[3]

v11 = pd.read_csv(IN_V11_TX_GENE_TAU)
gtex = pd.read_csv(IN_GTEX_TAU)
v11["gene_base"] = v11["gene"].str.split(".").str[0]
gtex["gene_base"] = gtex["gene"].str.split(".").str[0]

df = v11.merge(gtex, on="gene_base", suffixes=("_v11txgene", "_gtex"))
print(f"genes matched by base id: {len(df)} (v11_tx_gene={len(v11)}, gtex={len(gtex)})", file=sys.stderr)

fig, axes = plt.subplots(1, 2, figsize=(14, 7))

for ax, col, label in [
    (axes[0], "tau_unweighted", "tau_unweighted"),
    (axes[1], "tau_free_weighted", "tau_free_weighted (hierclust0.1, max-agg)"),
]:
    d = df.dropna(subset=[f"{col}_v11txgene", f"{col}_gtex"])
    x = d[f"{col}_gtex"].to_numpy()
    y = d[f"{col}_v11txgene"].to_numpy()
    r = np.corrcoef(x, y)[0, 1]
    mae = np.mean(np.abs(x - y))
    frac_close = np.mean(np.abs(x - y) < 0.05)
    ax.scatter(x, y, s=6, alpha=0.1, color="steelblue")
    lims = [min(x.min(), y.min()), max(x.max(), y.max())]
    ax.plot(lims, lims, color="black", linestyle="--", linewidth=1)
    ax.set_xlabel(f"gtex_protein_coding {label}")
    ax.set_ylabel(f"gtex_v11_tx_gene {label}")
    ax.set_title(f"{label}\nn={len(d)} r={r:.4f} MAE={mae:.4f} |diff|<0.05: {frac_close*100:.1f}%")
    print(f"{col}: n={len(d)} pearson_r={r:.6f} MAE={mae:.6f} frac_|diff|<0.05={frac_close:.4f}",
          file=sys.stderr)

fig.suptitle("gtex_v11_tx_gene vs. gtex_protein_coding: gene-level tau agreement")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)
print(f"-> {OUT_PNG}", file=sys.stderr)
