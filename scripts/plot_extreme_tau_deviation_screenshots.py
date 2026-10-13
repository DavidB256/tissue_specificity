"""For each of the 91 extreme-tau-deviation transcripts (see
extract_extreme_tau_deviation_gtf.py), render an IGV-style "screenshot" of
its gene's full model (every GENCODE transcript of that gene, not just
quantified ones), highlighting the extreme transcript. One PNG per
transcript, into its own subdirectory."""
import sys
import os
import re
from collections import defaultdict
import pandas as pd
from gene_model_render import render_gene_model

IN_GTF = sys.argv[1]
IN_SELECTED_TX_CSV = sys.argv[2]
OUT_DIR = sys.argv[3]

selected = pd.read_csv(IN_SELECTED_TX_CSV)
selected_genes = set(selected["gene_id"])

genes = {}  # gene_id -> gene dict
gene_transcripts = defaultdict(dict)  # gene_id -> {transcript_id: {strand, exons, cds}}

with open(IN_GTF) as f:
    for line in f:
        if line[0] == "#":
            continue
        fields = line.rstrip("\n").split("\t", 8)
        feature = fields[2]
        attrs = fields[8]
        gene_id = re.search(r'gene_id "([^"]+)"', attrs).group(1)
        if gene_id not in selected_genes:
            continue

        if feature == "gene":
            gene_name = re.search(r'gene_name "([^"]+)"', attrs).group(1)
            genes[gene_id] = {
                "gene_id": gene_id, "gene_name": gene_name, "chrom": fields[0],
                "start": int(fields[3]), "end": int(fields[4]), "strand": fields[6],
            }
        elif feature in ("transcript", "exon", "CDS"):
            transcript_id = re.search(r'transcript_id "([^"]+)"', attrs).group(1)
            tx = gene_transcripts[gene_id].setdefault(
                transcript_id, {"strand": fields[6], "exons": [], "cds": []})
            start, end = int(fields[3]), int(fields[4])
            if feature == "exon":
                tx["exons"].append((start, end))
            elif feature == "CDS":
                tx["cds"].append((start, end))

os.makedirs(OUT_DIR, exist_ok=True)

n_written = 0
for _, row in selected.iterrows():
    gene_id, tx_id = row["gene_id"], row["transcript_id"]
    gene = genes[gene_id]
    transcripts = {t: d for t, d in gene_transcripts[gene_id].items() if d["exons"]}

    subtitle = (f"extreme tx {tx_id}: isoform_fraction={row['isoform_fraction']:.3f}, "
                f"tx_tau_weighted={row['tx_tau_weighted']:.3f}, gene_tau_weighted={row['gene_tau_weighted']:.3f}, "
                f"arctanh_deviation={row['deviation_arctanh']:.3f}")
    out_png = os.path.join(OUT_DIR, f"{gene['gene_name']}_{tx_id}.png")
    render_gene_model(gene, transcripts, tx_id, out_png, subtitle=subtitle)
    n_written += 1

print(f"n_screenshots={n_written} n_genes={len(genes)} -> {OUT_DIR}", file=sys.stderr)
