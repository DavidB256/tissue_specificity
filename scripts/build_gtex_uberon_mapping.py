"""GTEx tissue (normalized, matches gtex_v11_tx_protein_coding's GCT column
names) -> UBERON id, from GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt's
SMTSD/SMUBRID columns. 3 tissues have no UBERON SMUBRID (different ontologies
entirely -- GTEx annotates them with EFO/CL cell-line/cell-type terms, not
anatomy) and are manually mapped to their anatomical source instead, per
explicit instruction:
  Cells_Cultured_fibroblasts (EFO:0002009, "fibroblast") -> UBERON:0002067
    dermis -- GTEx's standard protocol cultures fibroblasts from skin punch
    biopsies, specifically the dermis layer.
  Cells_EBV-transformed_lymphocytes (EFO:0000572, lymphoblastoid cell line)
    -> UBERON:0000178 blood -- the cells' anatomical source.
  Liver_Hepatocyte (CL:0000182, "hepatocyte") -> UBERON:0001114 liver --
    reuses the SAME id as the Liver/Liver_Mixed_Cell tissues (a hepatocyte's
    anatomical source is the liver).
12 of the 68 tissue labels share only 6 distinct UBERON ids (e.g. Stomach and
Stomach_Mixed_Cell both -> UBERON:0000945) -- kept as-is (both labels get the
same uberon_id); compute_uberon_tissue_distance.py's pairwise distance comes
out 0 for same-id pairs, which average-linkage clustering naturally resolves
into sibling leaves, rather than needing special-cased tree surgery."""
import sys
import pandas as pd

IN_SAMPLE_ATTRS_TXT = sys.argv[1]
IN_GCT = sys.argv[2]
OUT_CSV = sys.argv[3]

MANUAL_OVERRIDE = {
    "Cells_Cultured_fibroblasts": "UBERON:0002067",
    "Cells_EBV-transformed_lymphocytes": "UBERON:0000178",
    "Liver_Hepatocyte": "UBERON:0001114",
}


def normalize_tissue(s):
    """Same scheme as build_gtex_v11_sample_tissue_map.normalize_tissue."""
    s = s.replace(" - ", "_")
    s = s.replace("(", "").replace(")", "")
    s = s.replace(" ", "_")
    return s


attrs = pd.read_csv(IN_SAMPLE_ATTRS_TXT, sep="\t", usecols=["SMTSD", "SMUBRID"], low_memory=False).drop_duplicates()
attrs["tissue"] = attrs["SMTSD"].map(normalize_tissue)

gct_tissues = list(pd.read_csv(IN_GCT, sep="\t", skiprows=2, nrows=0).columns[2:])
out = pd.DataFrame({"tissue": gct_tissues}).merge(attrs[["tissue", "SMUBRID"]], on="tissue", how="left")
out = out.rename(columns={"SMUBRID": "uberon_id"})
out["source"] = "smubrid"

for tissue, uberon_id in MANUAL_OVERRIDE.items():
    mask = out["tissue"] == tissue
    if mask.any():
        out.loc[mask, "uberon_id"] = uberon_id
        out.loc[mask, "source"] = "manual_override"

unmapped = out[out["uberon_id"].isna() | ~out["uberon_id"].str.startswith("UBERON")]
if len(unmapped):
    print(f"WARNING: {len(unmapped)} tissues still not mapped to a UBERON id:\n{unmapped}", file=sys.stderr)

out.to_csv(OUT_CSV, index=False)
print(f"tissues={len(out)} smubrid={sum(out['source'] == 'smubrid')} "
      f"manual_override={sum(out['source'] == 'manual_override')} "
      f"distinct_uberon_ids={out['uberon_id'].nunique()} -> {OUT_CSV}", file=sys.stderr)
