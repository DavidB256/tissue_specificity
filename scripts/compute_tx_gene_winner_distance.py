"""Tissue x tissue Pearson correlation matrix of a txs_source's protein-coding
transcript x tissue TPM matrix, plus per-transcript winner-tissue analysis: the
tissue where each transcript is most highly expressed ("wins" it) vs. the tissue
where its gene (transcript TPMs summed per gene_id) wins, and the distance
(1 - Pearson r, from the corr matrix) between those two tissues. Ties broken to
the alphabetically-first tissue -- tissue columns are already sorted
alphabetically upstream (build_gtex_v11_tx_gct.py / build_gtex_txs_gct.py), so
plain argmax (numpy returns the first max on ties) suffices."""
import sys
import numpy as np
import pandas as pd

IN_GCT = sys.argv[1]
OUT_CORR_CSV = sys.argv[2]
OUT_DIST_CSV = sys.argv[3]

df = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(df.columns[2:])
assert tissues == sorted(tissues), "tissue columns must be alphabetically sorted for tie-breaking"

tx_mat = df[tissues].to_numpy()
corr = pd.DataFrame(np.corrcoef(tx_mat, rowvar=False), index=tissues, columns=tissues)
corr.to_csv(OUT_CORR_CSV)

gene_mat_df = df.groupby("Description")[tissues].sum()
gene_mat = gene_mat_df.to_numpy()

tissues_arr = np.array(tissues)
tx_winner = tissues_arr[np.argmax(tx_mat, axis=1)]
gene_winner_map = pd.Series(tissues_arr[np.argmax(gene_mat, axis=1)], index=gene_mat_df.index)
tx_gene_winner = gene_winner_map.loc[df["Description"]].to_numpy()

tissue_pos = {t: i for i, t in enumerate(tissues)}
corr_vals = corr.to_numpy()
dist = 1.0 - corr_vals[[tissue_pos[t] for t in tx_winner], [tissue_pos[t] for t in tx_gene_winner]]

out = pd.DataFrame({
    "transcript_id": df["Name"],
    "gene_id": df["Description"],
    "winning_tissue_tx": tx_winner,
    "winning_tissue_gene": tx_gene_winner,
    "dist": dist,
})
out.to_csv(OUT_DIST_CSV, index=False)

n_same = int((tx_winner == tx_gene_winner).sum())
print(f"tissues={len(tissues)} -> {OUT_CORR_CSV}", file=sys.stderr)
print(f"transcripts={len(out)} same_tissue_as_gene={n_same} ({n_same / len(out) * 100:.1f}%) -> {OUT_DIST_CSV}",
      file=sys.stderr)
