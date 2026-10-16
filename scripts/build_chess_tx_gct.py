"""Transcript x tissue median-TPM GCT for the CHESS bodysite expression table.
Raw format: transcript_id, data -- data is ";"-joined per-tissue blobs
"Tissue Name:min:..,25:..,50:..,75:..,max:..,mean:..,std:..,count:N". Per-tissue
median = the "50:" field (same role as the rest of this pipeline's per-tissue
median TPM). A transcript with no blob for a tissue is 0 there (never observed
in any sample of that tissue). Median<1 floored to 0 immediately after
extraction, same convention as build_gtex_v11_tx_gct.py. Joined to gene_id via
IN_FULL_META_CSV (parse_chess_gtf_meta.py), filtered to gene_type=="protein_coding"
(CHESS has no per-transcript biotype, so this is the gene-level filter only --
the closest available equivalent to gtex_v11_tx_protein_coding's stricter
transcript-biotype filter). Tissue names normalized with the same scheme used
everywhere else in this pipeline (build_gtex_v11_sample_tissue_map.normalize_tissue)."""
import sys
import re
import gzip
import pandas as pd

IN_RAW_TSV_GZ = sys.argv[1]
IN_FULL_META_CSV = sys.argv[2]
OUT_GCT = sys.argv[3]


def normalize_tissue(s):
    """Same scheme as build_gtex_v11_sample_tissue_map.normalize_tissue."""
    s = s.replace(" - ", "_")
    s = s.replace("(", "").replace(")", "")
    s = s.replace(" ", "_")
    return s

MEDIAN_RE = re.compile(r"(?:^|,)50:([0-9.eE+-]+)")

raw = pd.read_csv(IN_RAW_TSV_GZ, sep="\t")

tissues = set()
rows = []
for transcript_id, data in zip(raw["transcript_id"], raw["data"]):
    medians = {}
    for chunk in data.split(";"):
        name, rest = chunk.split(":", 1)
        m = MEDIAN_RE.search(rest)
        if m:
            medians[normalize_tissue(name)] = float(m.group(1))
    tissues.update(medians)
    rows.append((transcript_id, medians))

tissues = sorted(tissues)
print(f"transcripts={len(rows)} tissues={len(tissues)}", file=sys.stderr)

mat = pd.DataFrame([{t: medians.get(t, 0.0) for t in tissues} for _, medians in rows],
                    index=[tid for tid, _ in rows])
mat[mat < 1] = 0.0

meta = pd.read_csv(IN_FULL_META_CSV).set_index("transcript_id")
pc_tx = meta.index[meta["gene_type"] == "protein_coding"]
mat = mat.loc[mat.index.intersection(pc_tx)]

out = pd.DataFrame({"Name": mat.index, "Description": meta.loc[mat.index, "gene_id"].values})
out = pd.concat([out.reset_index(drop=True), mat.reset_index(drop=True)], axis=1)

with gzip.open(OUT_GCT, "wt") as f:
    f.write("#1.2\n")
    f.write(f"{len(out)}\t{len(tissues)}\n")
    out.to_csv(f, sep="\t", index=False)

print(f"protein_coding transcripts={len(out)} -> {OUT_GCT}", file=sys.stderr)
