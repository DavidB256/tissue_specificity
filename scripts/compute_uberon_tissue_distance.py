"""Pairwise UBERON-ontology distance between GTEx tissues (build_gtex_uberon_
mapping.py), for comparison (via Robinson-Foulds / cophenetic distance, see
compare_dendrogram_scenarios.py) against the transcriptome-correlation
dendrograms. UBERON is a DAG (most terms have >1 is_a parent), so rather than
force an artificial single-parent subtree, distance is computed directly on
the DAG: nxontology's batet similarity (an information-content-weighted,
multi-parent-aware edge-counting metric) between every pair of tissues'
UBERON ids, using ONLY is_a (rdfs:subClassOf) edges -- NOT part_of/develops_
from/etc, which reflect physical containment or lineage rather than taxonomic
generality. IC is computed from the full UBERON is_a graph's topology
(intrinsic IC, no external corpus needed), not just the subgraph touching our
68 tissues, so the informativeness stats aren't distorted by our small sample.
distance = 1 - batet, in [0, 1]. Two tissue labels that map to the SAME
UBERON id (e.g. Stomach / Stomach_Mixed_Cell) get distance exactly 0 --
average-linkage clustering naturally merges them as sibling leaves first,
no special-casing needed. Output: a tissue x tissue distance CSV (same shape
as tissue_correlation.csv elsewhere in this pipeline) and a linkage CSV
(average linkage, same format as tissue_correlation_linkage.csv)."""
import sys
import warnings
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform

IN_UBERON_OWL = sys.argv[1]
IN_MAPPING_CSV = sys.argv[2]
OUT_DISTANCE_CSV = sys.argv[3]
OUT_LINKAGE_CSV = sys.argv[4]

warnings.filterwarnings("ignore")
import pronto
import networkx as nx
from nxontology import NXOntology

mapping = pd.read_csv(IN_MAPPING_CSV)
tissues = list(mapping["tissue"])
uberon_of = dict(zip(mapping["tissue"], mapping["uberon_id"]))
needed_ids = set(mapping["uberon_id"])

print(f"loading {IN_UBERON_OWL} ...", file=sys.stderr)
onto = pronto.Ontology(IN_UBERON_OWL)
print(f"loaded: {len(onto)} terms", file=sys.stderr)

missing = [u for u in needed_ids if u not in onto]
if missing:
    raise ValueError(f"UBERON ids not found in ontology: {missing}")

g = nx.DiGraph()
for term in onto.terms():
    if term.obsolete:
        continue
    g.add_node(term.id, name=term.name)
    for parent in term.superclasses(distance=1, with_self=False):
        if parent.obsolete:
            continue
        g.add_edge(parent.id, term.id)
print(f"is_a graph: {g.number_of_nodes()} nodes, {g.number_of_edges()} edges, "
      f"{nx.number_weakly_connected_components(g)} components", file=sys.stderr)

nxo = NXOntology(g)
nxo.freeze()

n = len(tissues)
dist = np.zeros((n, n))
n_errors = 0
for i in range(n):
    for j in range(i + 1, n):
        id_i, id_j = uberon_of[tissues[i]], uberon_of[tissues[j]]
        if id_i == id_j:
            d = 0.0
        else:
            try:
                d = 1.0 - nxo.similarity(id_i, id_j).batet
            except Exception as e:
                print(f"WARNING: similarity({id_i}, {id_j}) failed ({e}), using distance=1", file=sys.stderr)
                d = 1.0
                n_errors += 1
        dist[i, j] = dist[j, i] = d
if n_errors:
    print(f"{n_errors} pairs failed similarity computation (likely different weakly-connected components)",
          file=sys.stderr)

pd.DataFrame(dist, index=tissues, columns=tissues).to_csv(OUT_DISTANCE_CSV)

Z = linkage(squareform(dist, checks=False), method="average")
pd.DataFrame(Z, columns=["cluster1", "cluster2", "distance", "n_obs"]).to_csv(OUT_LINKAGE_CSV, index=False)

print(f"tissues={n} -> {OUT_DISTANCE_CSV}, {OUT_LINKAGE_CSV}", file=sys.stderr)
