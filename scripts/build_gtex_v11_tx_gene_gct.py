"""Gene-level protein-coding TPM GCT for gtex_v11_tx: sum transcript TPM per
gene_id (Description) within each tissue column, from the already TPM<1->0
floored transcript-level GCT (build_gtex_v11_tx_gct.py floors immediately after
the per-tissue median). Floor is not reapplied here -- a gene's per-tissue sum
of floored transcript values is already either 0 (every transcript <1 in that
tissue) or >=1 (at least one transcript >=1 in that tissue), so an extra floor
would be a no-op."""
import sys
import gzip
import pandas as pd

IN_TX_GCT = sys.argv[1]
OUT_GENE_GCT = sys.argv[2]

df = pd.read_csv(IN_TX_GCT, sep="\t", skiprows=2)
tissues = list(df.columns[2:])

gene_mat = df.groupby("Description")[tissues].sum()
out = gene_mat.reset_index().rename(columns={"Description": "Name"})
out.insert(1, "Description", out["Name"])
out = out[["Name", "Description"] + tissues]

with gzip.open(OUT_GENE_GCT, "wt") as f:
    f.write("#1.2\n")
    f.write(f"{len(out)}\t{len(tissues)}\n")
    out.to_csv(f, sep="\t", index=False)

print(f"genes={len(out)} (from {len(df)} transcripts) tissues={len(tissues)} -> {OUT_GENE_GCT}", file=sys.stderr)
