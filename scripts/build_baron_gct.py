"""Transpose the Baron pancreas cell-type x gene TPM matrix
(build_baron_celltype_tpm.py, rows=cell type, cols=gene) into a genes x
cell-types GCT, so it can feed the existing tissue-correlation scripts
(compute_tissue_correlation.py / plot_tissue_correlation.py) unmodified.
Name/Description are both the gene symbol -- Baron has no Ensembl IDs."""
import sys
import gzip
import pandas as pd

IN_CSV = sys.argv[1]
OUT_GCT = sys.argv[2]

df = pd.read_csv(IN_CSV, index_col=0)  # cell_type x gene
mat = df.T  # gene x cell_type
cell_types = list(mat.columns)

out = pd.DataFrame({"Name": mat.index, "Description": mat.index})
for ct in cell_types:
    out[ct] = mat[ct].values

with gzip.open(OUT_GCT, "wt") as f:
    f.write("#1.2\n")
    f.write(f"{len(out)}\t{len(cell_types)}\n")
    out.to_csv(f, sep="\t", index=False)

print(f"genes={len(out)} cell_types={len(cell_types)}", file=sys.stderr)
