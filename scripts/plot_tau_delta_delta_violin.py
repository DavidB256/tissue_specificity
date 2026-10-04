"""Violin plots of DeltaDeltaTau := |delta_tau_w| - |delta_tau| (negative = free
weighting shrank the tau shift from collapsing tissue resolution, positive = it grew),
grouped by each gene's most-expressed tissue cluster in hierclust0.1 (dist_threshold=0.1
hierarchical clustering, see identify tissue_clusters.csv). Two rows: onebrain leg,
hierclust0.1 leg. Same cluster assignment (from the full 68-tissue GCT) is used for
both rows. TPM ties for a gene's max-expressing tissue are broken randomly. Each violin is
annotated with n and the p-value of a one-sample t-test that the mean is 0."""
import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import ttest_1samp
from gls_tau import floor_tpm
from plot_utils import display_label

IN_GCT = sys.argv[1]
IN_CLUSTERS = sys.argv[2]
IN_ONEBRAIN_FREE_TAU = sys.argv[3]
IN_ONEBRAIN_REF_TAU = sys.argv[4]
IN_SIMPLIFIED_FREE_TAU = sys.argv[5]
IN_SIMPLIFIED_REF_TAU = sys.argv[6]
OUT_PNG = sys.argv[7]

rng = np.random.default_rng(0)

SCENARIOS = [
    ("onebrain", IN_ONEBRAIN_FREE_TAU, IN_ONEBRAIN_REF_TAU),
    ("hierclust0.1", IN_SIMPLIFIED_FREE_TAU, IN_SIMPLIFIED_REF_TAU),
]


def load_delta_delta(free_tau_csv, ref_tau_csv):
    d = pd.read_csv(free_tau_csv)[["gene", "tau_unweighted", "tau_free_weighted"]]
    r = pd.read_csv(ref_tau_csv)[["gene", "tau_unweighted"]].rename(columns={"tau_unweighted": "tau_ref"})
    df = d.merge(r, on="gene").dropna()
    delta_tau = (df["tau_unweighted"] - df["tau_ref"]).abs()
    delta_tau_w = (df["tau_free_weighted"] - df["tau_ref"]).abs()
    return pd.DataFrame({"gene": df["gene"], "delta_delta_tau": delta_tau_w - delta_tau})


# per-gene most-expressed tissue -> hierclust0.1 cluster, ties broken randomly
gct = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(gct.columns[2:])
mat = floor_tpm(gct[tissues].values.astype(float))
is_max = mat == mat.max(axis=1, keepdims=True)
argmax_idx = np.array([rng.choice(np.flatnonzero(row)) for row in is_max])
argmax_tissue = np.array(tissues)[argmax_idx]

clusters = pd.read_csv(IN_CLUSTERS)
tissue_to_cluster = dict(zip(clusters["tissue"], clusters["cluster_id"]))
cluster_label = (clusters.assign(namelen=clusters["tissue"].str.len())
                  .sort_values(["namelen", "tissue"])
                  .groupby("cluster_id")["tissue"].first().to_dict())

gene_cluster = pd.DataFrame({
    "gene": gct["Name"].values,
    "cluster": [cluster_label[tissue_to_cluster[t]] for t in argmax_tissue],
})

cluster_order = clusters.sort_values("cluster_id")["cluster_id"].map(cluster_label).drop_duplicates().tolist()

fig, axes = plt.subplots(2, 1, figsize=(max(12, 0.4 * len(cluster_order)), 12))
for ax, (label, free_tau_csv, ref_tau_csv) in zip(axes, SCENARIOS):
    dd = load_delta_delta(free_tau_csv, ref_tau_csv).merge(gene_cluster, on="gene")
    counts = dd["cluster"].value_counts()
    order = [c for c in cluster_order if c in counts.index]

    sns.violinplot(data=dd, x="cluster", y="delta_delta_tau", order=order, ax=ax,
                    color="tab:blue", cut=0, inner="quartile")
    ax.axhline(0, color="gray", lw=1, ls="--")
    ax.set_xlabel("tissue cluster (labeled by shortest tissue name, hierclust0.1)")
    ax.set_ylabel(r"$\Delta\Delta\tau = |\Delta\tau_w| - |\Delta\tau|$")
    ax.set_title(f"full vs {display_label(label)}, n={len(dd)} genes", fontsize=10, pad=14)
    ax.tick_params(axis="x", rotation=30, labelsize=7)
    for tick in ax.get_xticklabels():
        tick.set_ha("right")

    ymin, ymax = dd["delta_delta_tau"].min(), dd["delta_delta_tau"].max()
    pad = 0.06 * (ymax - ymin)
    ax.set_ylim(ymin - pad, ymax + 3 * pad)
    for i, c in enumerate(order):
        vals = dd.loc[dd["cluster"] == c, "delta_delta_tau"]
        pval = ttest_1samp(vals, 0).pvalue
        ax.text(i, vals.max() + pad * 0.3, f"n={counts[c]}\np={pval:.1e}",
                ha="center", va="bottom", fontsize=6, rotation=30)

fig.suptitle(r"$\Delta\Delta\tau$ by gene's most-expressed tissue cluster (hierclust0.1)")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)
