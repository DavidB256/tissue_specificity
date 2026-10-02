"""Gene-level tau_unweighted for gtex_txs: aggregate the transcript x tissue
protein-coding TPM matrix (base_gct) up to gene x tissue by summing transcript TPM
per gene_id (Description col) within each tissue -- the standard way to recover gene
expression from transcript-level quantification. No TPM<1 flooring (see
compute_tau_unweighted.py). Also reports median_log_tpm: each gene's median (across
the 9 tissues) of log10(TPM+1) -- a single overall-expression-level scalar, used for
coloring points in plot_gene_vs_transcript_tau_variance.py."""
import sys
import numpy as np
import pandas as pd
from gls_tau import weighted_tau

IN_GCT = sys.argv[1]
OUT_CSV = sys.argv[2]

df = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(df.columns[2:])

gene_mat = df.groupby("Description")[tissues].sum()
mat = gene_mat.to_numpy()

median_log_tpm = np.median(np.log10(mat + 1.0), axis=1)

out = pd.DataFrame({
    "gene_id": gene_mat.index,
    "tau_unweighted": weighted_tau(mat, np.ones(len(tissues))),
    "median_log_tpm": median_log_tpm,
})
out.to_csv(OUT_CSV, index=False)

print(f"genes={len(out)} (from {len(df)} transcripts) tissues={len(tissues)} -> {OUT_CSV}", file=sys.stderr)
