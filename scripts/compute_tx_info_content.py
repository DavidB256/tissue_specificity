"""Per-gene "information conveyed about which transcript is expressed": Shannon
entropy H of the gene's isoform-usage distribution p_i = tx_i's total TPM
(summed across tissues) / gene's total TPM, restricted to transcripts with
nonzero total TPM (an all-zero transcript has no observed tissue profile to
define a KL term against, and contributes nothing to H either way). H=0 (all
info in one isoform) for a single-expressed-transcript gene; H increases
toward log2(n) as usage flattens across more, more-evenly-used isoforms.
H_norm = 1 - H/log2(n) (defined as 1 for n=1) bounds this to [0,1] regardless
of n. effective_n = 2**H is the "effective number of isoforms in play"
(perplexity).

weighted_kl: for each transcript, KL(tx_tissue_pmf || gene_tissue_pmf) in bits
-- always finite here, since the gene's per-tissue TPM is the SUM of its
transcripts', so gene_pmf[t]==0 implies every transcript's pmf is 0 there too
(never log(positive/0)). weighted_kl = sum_i p_i * KL_i, the TPM-share-
weighted average divergence of each isoform's tissue-usage pattern from the
gene's overall pattern. Equivalently: the mutual information between "which
transcript" and "which tissue" under the model tx ~ p_i, tissue ~
tx's own tissue pmf.

Writes two outputs: OUT_PICKLE, a dict gene_id -> {transcript_id: (tx_total_tpm,
KL_i)} for every expressed transcript, and OUT_CSV with one row per gene:
gene_id, n_tx (expressed transcripts used), H, H_norm, effective_n,
weighted_kl. Genes with zero expressed transcripts (every isoform floored to
0 TPM in every tissue) are dropped -- no distribution to measure."""
import sys
import pickle
import numpy as np
import pandas as pd

IN_GCT = sys.argv[1]
OUT_PICKLE = sys.argv[2]
OUT_CSV = sys.argv[3]

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

out_dict = {}
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

    p = tx_total[idx] / gene_total
    n = len(idx)
    H = -float(np.sum(p * np.log2(p)))
    H_norm = 1.0 if n == 1 else 1.0 - H / np.log2(n)
    effective_n = 2.0 ** H

    gene_entry = {}
    weighted_kl = 0.0
    for j, i in enumerate(idx):
        mask = tx_pmf[i] > 0
        ratio = np.divide(tx_pmf[i], gene_pmf, out=np.ones_like(tx_pmf[i]), where=mask)
        kl = float(np.sum(np.where(mask, tx_pmf[i] * np.log2(ratio), 0.0)))
        gene_entry[tx_ids[i]] = (float(tx_total[i]), kl)
        weighted_kl += p[j] * kl
    out_dict[gene_id] = gene_entry

    rows.append({
        "gene_id": gene_id, "n_tx": n, "H": H, "H_norm": H_norm,
        "effective_n": effective_n, "weighted_kl": weighted_kl,
    })

with open(OUT_PICKLE, "wb") as f:
    pickle.dump(out_dict, f)

out = pd.DataFrame(rows)
out.to_csv(OUT_CSV, index=False)

print(f"transcripts: {len(df)} total, {n_dropped_tx} dropped (zero TPM in every tissue)", file=sys.stderr)
print(f"genes: {out['gene_id'].nunique() if len(out) else 0} written, "
      f"{n_genes_no_expr} dropped (no expressed transcripts) -> {OUT_CSV}", file=sys.stderr)
print(f"-> {OUT_PICKLE}", file=sys.stderr)
