"""Concise 3-panel figure for compute_uberon_mantel.py's output: one scatter
per source, UBERON anatomical distance (1-batet) vs. transcriptomic
cophenetic distance over all common-tissue pairs, annotated with the Mantel
r and permutation p-value. A flat/scattered cloud with r near 0 means
anatomical (UBERON) relatedness does not track transcriptomic similarity for
that source."""
import sys
import pickle
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IN_RESULTS_CSV = sys.argv[1]
OUT_PNG = sys.argv[2]

SOURCE_LABEL = {
    "gtex_protein_coding": "GENCODE (gene)",
    "gtex_v11_tx_protein_coding": "GENCODE (tx)",
    "chess_tx_protein_coding": "CHESS (tx)",
}

results = pd.read_csv(IN_RESULTS_CSV).set_index("source")
with open(IN_RESULTS_CSV.replace(".csv", "_vectors.pkl"), "rb") as f:
    mantel_data = pickle.load(f)

sources = list(results.index)
fig, axes = plt.subplots(1, len(sources), figsize=(6 * len(sources), 5.5), sharey=True)
for ax, name in zip(axes, sources):
    d = mantel_data[name]
    r = results.loc[name, "mantel_r"]
    p = results.loc[name, "p_value"]
    n = int(results.loc[name, "n_common_tissues"])
    ax.scatter(d["uberon_dist"], d["coph_dist"], s=6, alpha=0.25, color="tab:blue")
    ax.set_xlabel("UBERON distance (1 - batet)")
    ax.set_title(SOURCE_LABEL.get(name, name), fontsize=11)
    ax.text(0.03, 0.97, f"mantel r={r:.3f}\np={p:.4f}\nn tissues={n}",
            transform=ax.transAxes, va="top", ha="left", fontsize=9,
            bbox=dict(facecolor="white", alpha=0.7, edgecolor="none"))
axes[0].set_ylabel("transcriptomic cophenetic distance")

fig.suptitle("UBERON anatomical distance vs. transcriptomic cophenetic distance (Mantel test)")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)

print(f"-> {OUT_PNG}", file=sys.stderr)
