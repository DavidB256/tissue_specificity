"""Free-weights analog of tau_weighted at a reduced (clustered) resolution, with no
GLS involved anywhere: collapse the per-tissue free_weights vector (fit on the default
resolution against this same reference, see optimize_free_weights.py) into one weight
per cluster -- sum of its member tissues' weights, matching build_gct_clustered.py's
mean-TPM collapse (so a cluster's total weight budget stays comparable to a single
independent tissue's) -- then compute weighted_tau on the already-collapsed GCT using
those per-cluster weights. Output has the same schema as optimize_free_weights.py's
tau output (gene, tau_unweighted, tau_free_weighted) so both resolutions can be
compared like-for-like."""
import sys
import numpy as np
import pandas as pd
from gls_tau import weighted_tau, floor_tpm

IN_GCT = sys.argv[1]
IN_WEIGHTS = sys.argv[2]
IN_CLUSTERS = sys.argv[3]
OUT_CSV = sys.argv[4]

gct = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(gct.columns[2:])
genes = gct["Name"].values
mat = floor_tpm(gct[tissues].values.astype(float))

weights = pd.read_csv(IN_WEIGHTS).set_index("tissue")["weight"]
clusters = pd.read_csv(IN_CLUSTERS)

collapsed_weight = {}
for cid, grp in clusters.groupby("cluster_id"):
    members = grp["tissue"].tolist()
    name = members[0] if len(members) == 1 else f"cluster{cid}"
    collapsed_weight[name] = weights[members].sum()

w = np.array([collapsed_weight[t] for t in tissues])
print(f"collapsed weights: {dict(zip(tissues, w))}", file=sys.stderr)

out = pd.DataFrame({
    "gene": genes,
    "tau_unweighted": weighted_tau(mat, np.ones(len(tissues))),
    "tau_free_weighted": weighted_tau(mat, w),
})
out.to_csv(OUT_CSV, index=False)
print(f"genes={len(out)} tissues={len(tissues)}", file=sys.stderr)
