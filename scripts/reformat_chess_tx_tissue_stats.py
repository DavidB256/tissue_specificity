"""Reformat the CHESS bodysite raw table (one row per transcript, one "data" col
of ";"-joined per-tissue blobs "Tissue:min:..,25:..,50:..,75:..,max:..,mean:..,
std:..,count:N") into a wide CSV: one row per transcript, one column per
(tissue, statistic) pair -- "{tissue}_{stat}" for stat in
min/25/50/75/max/mean/std/count. A transcript with no blob for a tissue has
NaN in all of that tissue's columns (never observed in any sample there, vs.
build_chess_tx_gct.py's median pipeline where an unobserved tissue is TPM=0 --
here we keep the distinction, since NaN vs 0 matters for a parametric bootstrap
fit). Tissue names normalized with the same scheme as the rest of this
pipeline (build_gtex_v11_sample_tissue_map.normalize_tissue)."""
import sys
import pandas as pd

IN_RAW_TSV_GZ = sys.argv[1]
OUT_CSV = sys.argv[2]

STATS = ["min", "25", "50", "75", "max", "mean", "std", "count"]


def normalize_tissue(s):
    """Same scheme as build_gtex_v11_sample_tissue_map.normalize_tissue."""
    s = s.replace(" - ", "_")
    s = s.replace("(", "").replace(")", "")
    s = s.replace(" ", "_")
    return s


raw = pd.read_csv(IN_RAW_TSV_GZ, sep="\t")

rows = []
tissues = set()
for transcript_id, data in zip(raw["transcript_id"], raw["data"]):
    row = {"transcript_id": transcript_id}
    for chunk in data.split(";"):
        name, rest = chunk.split(":", 1)
        tissue = normalize_tissue(name)
        tissues.add(tissue)
        for kv in rest.split(","):
            key, val = kv.split(":", 1)
            row[f"{tissue}_{key}"] = float(val)
    rows.append(row)

tissues = sorted(tissues)
columns = ["transcript_id"] + [f"{t}_{stat}" for t in tissues for stat in STATS]

out = pd.DataFrame(rows, columns=columns)
out.to_csv(OUT_CSV, index=False)

print(f"transcripts={len(out)} tissues={len(tissues)} cols={len(out.columns)} -> {OUT_CSV}", file=sys.stderr)
