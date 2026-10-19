"""Per-gene mean Wasserstein-1 distance between each expressed transcript's
tissue PMF and the gene's own (summed-isoform) tissue PMF -- same PMF/
Wasserstein convention as find_high_tau_isoform_pairs_wasserstein.py
(tissues ordered by fixed GCT column index, shared support 0..n_tissues-1).
Unweighted (simple) mean across a gene's expressed transcripts -- a
single-isoform gene's tx PMF equals its gene PMF exactly, so it gets
mean_wasserstein=0 trivially. Genes with zero expressed transcripts (every
isoform floored to 0 TPM in every tissue) are dropped."""
import sys
import numpy as np
import pandas as pd
from scipy.stats import wasserstein_distance

IN_GCT = sys.argv[1]
OUT_CSV = sys.argv[2]

df = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(df.columns[2:])
mat = df[tissues].to_numpy(dtype=float)
tx_ids = df["Name"].to_numpy()
gene_ids = df["Description"].to_numpy()

tx_total = mat.sum(axis=1)
keep = tx_total > 0
n_dropped_tx = int((~keep).sum())
mat, tx_ids, gene_ids, tx_total = mat[keep], tx_ids[keep], gene_ids[keep], tx_total[keep]
tx_pmf = mat / tx_total[:, None]

positions = np.arange(len(tissues), dtype=float)

rows = []
n_genes_no_expr = 0
for gene_id in pd.unique(gene_ids):
    idx = np.where(gene_ids == gene_id)[0]
    gene_profile = mat[idx].sum(axis=0)
    gene_total = gene_profile.sum()
    if gene_total <= 0:
        n_genes_no_expr += 1
        continue
    gene_pmf = gene_profile / gene_total

    dists = [wasserstein_distance(positions, positions, u_weights=tx_pmf[i], v_weights=gene_pmf)
             for i in idx]
    rows.append({"gene_id": gene_id, "n_tx": len(idx), "mean_wasserstein_tx_to_gene": float(np.mean(dists))})

out = pd.DataFrame(rows)
out.to_csv(OUT_CSV, index=False)

print(f"transcripts: {len(df)} total, {n_dropped_tx} dropped (zero TPM in every tissue)", file=sys.stderr)
print(f"genes: {len(out)} written, {n_genes_no_expr} dropped (no expressed transcripts) -> {OUT_CSV}",
      file=sys.stderr)
