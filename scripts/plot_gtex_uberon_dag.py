"""DAG of all 68 GTEx tissues' UBERON terms, together with every is_a ancestor
needed to connect them (since most tissue-pairs share only a distant common
ancestor, not a direct edge) -- same is_a-only, obsolete-filtered graph as
compute_uberon_tissue_distance.py, just rendered directly rather than reduced
to a pairwise-distance tree. Leaf nodes (actual GTEx tissues, from
build_gtex_uberon_mapping.py) are labeled with their GTEx tissue name(s) --
multiple names when several GTEx labels share one UBERON id, e.g. Stomach /
Stomach_Mixed_Cell -- plus the UBERON id; ancestor-only nodes are labeled with
just the UBERON term name and id."""
import sys
import warnings
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

warnings.filterwarnings("ignore")
import pronto

IN_UBERON_OWL = sys.argv[1]
IN_MAPPING_CSV = sys.argv[2]
OUT_PNG = sys.argv[3]

LEAF_COLOR = "#4C72B0"
ANCESTOR_COLOR = "#DDDDDD"
LEAF_FONTSIZE = 6.5
ANCESTOR_FONTSIZE = 5

mapping = pd.read_csv(IN_MAPPING_CSV)
tissues_of = mapping.groupby("uberon_id")["tissue"].apply(list).to_dict()
leaf_ids = set(tissues_of)

print(f"loading {IN_UBERON_OWL} ...", file=sys.stderr)
onto = pronto.Ontology(IN_UBERON_OWL)

g = nx.DiGraph()
for term in onto.terms():
    if term.obsolete:
        continue
    for parent in term.superclasses(distance=1, with_self=False):
        if parent.obsolete:
            continue
        g.add_edge(parent.id, term.id)

needed = set(leaf_ids)
for lid in leaf_ids:
    needed |= nx.ancestors(g, lid)
sub = g.subgraph(needed).copy()
print(f"DAG: {sub.number_of_nodes()} nodes ({len(leaf_ids)} GTEx tissue terms, "
      f"{sub.number_of_nodes() - len(leaf_ids)} ancestors-only), {sub.number_of_edges()} edges",
      file=sys.stderr)

safe_of = {n: n.replace(":", "_") for n in sub.nodes}
safe_graph = nx.relabel_nodes(sub, safe_of, copy=True)
for n in safe_graph.nodes:
    safe_graph.nodes[n]["label"] = ""  # suppress graphviz's own text sizing from affecting layout
safe_graph.graph["graph"] = {"nodesep": "0.7", "ranksep": "0.9"}
pos_safe = nx.nx_pydot.graphviz_layout(
    safe_graph, prog="dot",
    root=None,
)
inv = {v: k for k, v in safe_of.items()}
pos = {inv[k]: v for k, v in pos_safe.items()}

xs = [p[0] for p in pos.values()]
ys = [p[1] for p in pos.values()]
width_in = (max(xs) - min(xs)) / 72 + 4
height_in = (max(ys) - min(ys)) / 72 + 3
print(f"figure size: {width_in:.0f} x {height_in:.0f} in", file=sys.stderr)

fig, ax = plt.subplots(figsize=(width_in, height_in))

nx.draw_networkx_edges(sub, pos, ax=ax, arrows=True, arrowsize=8, edge_color="#999999", width=0.6, node_size=0)

leaf_nodes = [n for n in sub.nodes if n in leaf_ids]
ancestor_nodes = [n for n in sub.nodes if n not in leaf_ids]
nx.draw_networkx_nodes(sub, pos, nodelist=ancestor_nodes, ax=ax, node_color=ANCESTOR_COLOR,
                        node_size=60, edgecolors="#888888", linewidths=0.4)
nx.draw_networkx_nodes(sub, pos, nodelist=leaf_nodes, ax=ax, node_color=LEAF_COLOR,
                        node_size=90, edgecolors="black", linewidths=0.6)

for n in ancestor_nodes:
    x, y = pos[n]
    ax.text(x, y - 9, f"{onto[n].name}\n{n}", ha="center", va="top", fontsize=ANCESTOR_FONTSIZE, color="#555555")

for n in leaf_nodes:
    x, y = pos[n]
    names = ",\n".join(sorted(tissues_of[n]))
    ax.text(x, y - 9, f"{names}\n{n}", ha="center", va="top", fontsize=LEAF_FONTSIZE,
            fontweight="bold", color="black")

legend_elems = [
    Line2D([0], [0], marker="o", color="none", markerfacecolor=LEAF_COLOR, markeredgecolor="black",
           markersize=9, label="GTEx tissue (68 total, some share a UBERON id)"),
    Line2D([0], [0], marker="o", color="none", markerfacecolor=ANCESTOR_COLOR, markeredgecolor="#888888",
           markersize=7, label="UBERON ancestor (not a GTEx tissue)"),
]
ax.legend(handles=legend_elems, loc="upper right", fontsize=10)

ax.set_title("UBERON is_a DAG spanning all 68 GTEx tissues", fontsize=14)
ax.set_xlim(min(xs) - 60, max(xs) + 60)
ax.set_ylim(min(ys) - 40, max(ys) + 40)
ax.axis("off")

fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)

print(f"-> {OUT_PNG}", file=sys.stderr)
