"""Cut the tissue hierarchical-clustering dendrogram (from compute_tissue_correlation.py)
at t=N_BLOCKS (first N_BLOCKS-1 splits) into tissue blocks. Used for cluster-level
resampling in compute_tau_bootstrap_variance.py; also recorded here as an intermediate
data file so the block membership is inspectable independent of that script."""
import sys
import pandas as pd
from scipy.cluster.hierarchy import fcluster

IN_CORR = sys.argv[1]
IN_LINKAGE = sys.argv[2]
N_BLOCKS = int(sys.argv[3])
OUT_CSV = sys.argv[4]

tissues = list(pd.read_csv(IN_CORR, index_col=0).index)
Z = pd.read_csv(IN_LINKAGE).values
labels = fcluster(Z, t=N_BLOCKS, criterion="maxclust")

out = pd.DataFrame({"tissue": tissues, "block": labels})
out.to_csv(OUT_CSV, index=False)

for b in sorted(out["block"].unique()):
    members = out.loc[out["block"] == b, "tissue"].tolist()
    print(f"block {b} ({len(members)}): {members}", file=sys.stderr)
