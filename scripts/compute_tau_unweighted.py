"""Plain unweighted tau (Yanai et al. 2005) only -- no diag_reg_inv_corr/margin/free
weighting. Same output schema (gene,tau_unweighted) as the other sources'
weighted_tau.csv, minus the weighted column, for gtex_txs (no leg/reference exists to
weight against -- see gtex-txs-pipeline plan). No TPM<1 flooring: unlike the gene-
level GTEx/Baron panels, every row here is a CHESS-annotated transcript, not a
candidate-call that might be noise, so a low-but-real TPM isn't treated as zero."""
import sys
import numpy as np
import pandas as pd
from gls_tau import weighted_tau

IN_GCT = sys.argv[1]
OUT_CSV = sys.argv[2]

df = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(df.columns[2:])
mat = df[tissues].values.astype(float)

out = pd.DataFrame({
    "gene": df["Name"].values,
    "tau_unweighted": weighted_tau(mat, np.ones(len(tissues))),
})
out.to_csv(OUT_CSV, index=False)

print(f"transcripts={len(out)} tissues={len(tissues)} -> {OUT_CSV}", file=sys.stderr)
