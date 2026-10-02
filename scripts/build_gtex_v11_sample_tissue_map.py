"""SAMPID -> tissue map for GTEx v11, from GTEx_Analysis_v11_Annotations_
SampleAttributesDS.txt's SAMPID/SMTSD columns. Tissue names are SMTSD normalized to
the same scheme GTEx itself uses for the gene_median_tpm GCT's column headers (" - "
-> "_", parens dropped but their contents kept, remaining spaces -> "_") -- verified
to reproduce all 68 of that GCT's tissue names exactly, so results here are directly
comparable/joinable to the rest of this pipeline's GTEx tissue naming."""
import re
import sys
import pandas as pd

IN_SAMPLE_ATTRS_TXT = sys.argv[1]
OUT_CSV = sys.argv[2]


def normalize_tissue(s):
    s = s.replace(" - ", "_")
    s = s.replace("(", "").replace(")", "")
    s = s.replace(" ", "_")
    return s


attrs = pd.read_csv(IN_SAMPLE_ATTRS_TXT, sep="\t", usecols=["SAMPID", "SMTSD"], low_memory=False)
attrs = attrs.dropna(subset=["SMTSD"]).drop_duplicates("SAMPID")
out = pd.DataFrame({
    "sample_id": attrs["SAMPID"].values,
    "tissue": attrs["SMTSD"].map(normalize_tissue).values,
})
out.to_csv(OUT_CSV, index=False)

print(f"samples={len(out)} tissues={out['tissue'].nunique()} -> {OUT_CSV}", file=sys.stderr)
