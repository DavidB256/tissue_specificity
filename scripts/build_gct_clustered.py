"""Collapse tissues into one column per hierarchical cluster (raw TPM, pre floor_tpm),
using cluster assignments from compute_tau_cluster_weighted.py's tissue_clusters.csv
(average-linkage clustering, distance=1-corr, cut at a fixed DIST_THRESHOLD -- e.g.
cluster_weights/full/tissue_clusters.csv for threshold=0.1). Singleton clusters
(cluster_size==1) keep their original tissue name; multi-tissue clusters become
"cluster{cluster_id}", collapsed via AGG_METHOD (mean, or max to keep the peak value
instead of diluting toward the cross-tissue average). Output is a standard 2-line-header
GCT (matches the input format)."""
import sys
import gzip
import pandas as pd

IN_GCT = sys.argv[1]
IN_CLUSTERS_CSV = sys.argv[2]
OUT_GCT = sys.argv[3]
AGG_METHOD = sys.argv[4] if len(sys.argv) > 4 else "mean"
assert AGG_METHOD in ("mean", "max"), f"AGG_METHOD must be mean or max, got {AGG_METHOD!r}"

df = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
clusters = pd.read_csv(IN_CLUSTERS_CSV)

out = df[["Name", "Description"]].copy()
new_tissues = []
for cid, grp in clusters.groupby("cluster_id"):
    members = grp["tissue"].tolist()
    if len(members) == 1:
        name = members[0]
        out[name] = df[name]
    else:
        name = f"cluster{cid}"
        out[name] = getattr(df[members], AGG_METHOD)(axis=1)
        print(f"{name} <- {members}", file=sys.stderr)
    new_tissues.append(name)

new_tissues = sorted(new_tissues)
out = out[["Name", "Description"] + new_tissues]

print(f"tissues {len(clusters)} -> {len(new_tissues)}", file=sys.stderr)

with gzip.open(OUT_GCT, "wt") as f:
    f.write("#1.2\n")
    f.write(f"{len(out)}\t{len(new_tissues)}\n")
    out.to_csv(f, sep="\t", index=False)
