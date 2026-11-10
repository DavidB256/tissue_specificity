"""Genes with >MIN_HITS hit tx (from list_high_dev_isoforms.py): gene symbol, n_hits, gene tau_w
(weights applied to gene-summed profile), top tissue (argmax gene-summed TPM).
Usage: HITS_CSV GCT WEIGHTS_CSV GTF OUT_CSV [MIN_HITS=2]"""
import sys
import re
import pandas as pd
from tau_utils import weighted_tau

hits_csv, gct_path, w_csv, gtf, out = sys.argv[1:6]
MIN_HITS = int(sys.argv[6]) if len(sys.argv) > 6 else 2

c = pd.read_csv(hits_csv)["gene_id"].value_counts()
c = c[c > MIN_HITS]
gct = pd.read_csv(gct_path, sep="\t", skiprows=2)
t = list(gct.columns[2:])
w = pd.read_csv(w_csv).set_index("tissue").loc[t]["weight"].to_numpy()
gm = gct.groupby("Description")[t].sum().loc[c.index]

names = {}
with open(gtf) as f:
    for line in f:
        if line[0] != "#" and "\tgene\t" in line:
            m = re.search(r'gene_id "([^"]+)".*gene_name "([^"]+)"', line)
            if m and m[1] in c.index:
                names[m[1]] = m[2]

res = pd.DataFrame({
    "gene_name": [names.get(g, "") for g in c.index],
    "gene_id": c.index,
    "n_hits": c.to_numpy(),
    "gene_tau_w": weighted_tau(gm.to_numpy(), w),
    "top_tissue": gm.idxmax(axis=1).to_numpy(),
})
res.to_csv(out, index=False)
print(f"n_genes={len(res)} -> {out}", file=sys.stderr)
