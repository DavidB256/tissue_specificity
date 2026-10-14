"""For each gene with >=1 extreme-tau-deviation transcript (see
extract_extreme_tau_deviation_gtf.py), render an IGV-style "screenshot" of
its full gene model (every GENCODE transcript of that gene, not just
quantified ones), highlighting the extreme transcript(s). Transcripts
ordered top-to-bottom by total TPM (summed across tissues); unquantified
transcripts (0 TPM) sort to the bottom. Grey vertical lines mark the TSS and
5' splice site of each extreme transcript. One PNG per gene (not per tx)."""
import sys
import os
import re
from collections import defaultdict
import pandas as pd
from gene_model_render import render_gene_model

IN_GTF = sys.argv[1]
IN_GCT = sys.argv[2]
IN_SELECTED_TX_CSV = sys.argv[3]
OUT_DIR = sys.argv[4]

selected = pd.read_csv(IN_SELECTED_TX_CSV)
selected_genes = set(selected["gene_id"])
gene_to_highlight_tx = selected.groupby("gene_id")["transcript_id"].apply(set).to_dict()

gct = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(gct.columns[2:])
tx_tpm = dict(zip(gct["Name"], gct[tissues].sum(axis=1)))

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


def tss_and_five_prime_ss(tx):
    exons_sorted = sorted(tx["exons"])
    if tx["strand"] == "+":
        first = exons_sorted[0]
        tss = first[0]
        five_prime_ss = first[1] if len(exons_sorted) > 1 else None
    else:
        first = exons_sorted[-1]
        tss = first[1]
        five_prime_ss = first[0] if len(exons_sorted) > 1 else None
    return tss, five_prime_ss


os.makedirs(OUT_DIR, exist_ok=True)

n_written = 0
for gene_id, highlight_tx_ids in gene_to_highlight_tx.items():
    gene = genes[gene_id]
    transcripts = {t: d for t, d in gene_transcripts[gene_id].items() if d["exons"]}
    order = sorted(transcripts, key=lambda t: (-tx_tpm.get(t, 0.0), t))

    vlines = []
    for tx_id in sorted(highlight_tx_ids):
        tss, five_prime_ss = tss_and_five_prime_ss(transcripts[tx_id])
        vlines.append((tss, f"TSS ({tx_id})"))
        if five_prime_ss is not None:
            vlines.append((five_prime_ss, f"5'SS ({tx_id})"))

    rows = selected[selected["gene_id"] == gene_id]
    subtitle_parts = [
        f"{r.transcript_id}: isoform_fraction={r.isoform_fraction:.3f}, tx_tau_weighted={r.tx_tau_weighted:.3f}, "
        f"arctanh_deviation={r.deviation_arctanh:.3f}"
        for r in rows.itertuples()
    ]
    subtitle = "extreme tx -- " + " | ".join(subtitle_parts)

    out_png = os.path.join(OUT_DIR, f"{gene['gene_name']}_{gene_id}.png")
    render_gene_model(gene, transcripts, order, highlight_tx_ids, vlines, out_png, subtitle=subtitle)
    n_written += 1

print(f"n_screenshots={n_written} n_genes={len(genes)} n_selected_tx={len(selected)} -> {OUT_DIR}", file=sys.stderr)
