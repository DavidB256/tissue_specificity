"""Shared utilities for comparing hierarchical-clustering dendrograms (any
average-linkage scipy linkage Z + its leaf-label list): bipartitions
(Robinson-Foulds), scenario construction, and the full pairwise cophenetic +
RF + ARI comparison engine. Used by compare_dendrogram_scenarios.py and
compare_gene_vs_tx_dendrograms.py -- see the former's docstring for what each
metric means and why (cophenetic = height-sensitive, RF = height-blind
topology-only, both restricted to tissues common to the pair being compared,
since different scenarios can have different leaf sets)."""
from itertools import combinations
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, cophenet, to_tree
from scipy.spatial.distance import squareform
from scipy.stats import pearsonr
from sklearn.metrics import adjusted_rand_score


def get_bipartitions(Z, tissues):
    """Non-trivial bipartitions (frozenset of tissue names on the smaller side)
    induced by every internal node of a linkage tree."""
    tree = to_tree(Z)
    n = len(tissues)
    splits = set()

    def recurse(node):
        if node.is_leaf():
            return frozenset([tissues[node.id]])
        combined = recurse(node.left) | recurse(node.right)
        if 2 <= len(combined) <= n - 2:
            splits.add(combined)
        return combined

    recurse(tree)
    return splits


def canonical_split(s, universe):
    """A split is only meaningful up to complementation (S vs. universe-S
    represent the same partition) -- canonicalize to whichever side is
    smaller (ties broken by sorted membership) so the same split from two
    different trees compares equal."""
    complement = universe - s
    return min(frozenset(s), frozenset(complement), key=lambda x: (len(x), sorted(x)))


def restricted_splits(splits, common_set):
    out = set()
    for s in splits:
        r = s & common_set
        if 2 <= len(r) <= len(common_set) - 2:
            out.add(canonical_split(r, common_set))
    return out


def make_scenario(Z, tissues, dist_threshold=0.1):
    cluster_id = fcluster(Z, t=dist_threshold, criterion="distance")
    coph = squareform(cophenet(Z))
    return {
        "tissues": tissues, "cophenetic": coph,
        "cluster_of": dict(zip(tissues, cluster_id)),
        "splits": get_bipartitions(Z, tissues),
        "n_clusters": len(set(cluster_id)),
    }


def compare_all_scenarios(scenarios):
    """scenarios: dict name -> scenario dict (from make_scenario). Returns
    (dist_matrix_df, rf_matrix_df, long_df). A pair with 0 common leaves
    (e.g. a gene-centric source vs. Baron's cell-type names) gets NaN in
    both matrices -- genuinely no shared leaves to compare, not an error."""
    names = list(scenarios.keys())
    dist_matrix = pd.DataFrame(np.zeros((len(names), len(names))), index=names, columns=names)
    rf_matrix = pd.DataFrame(np.zeros((len(names), len(names))), index=names, columns=names)
    long_rows = []
    for a, b in combinations(names, 2):
        sa, sb = scenarios[a], scenarios[b]
        common = sorted(set(sa["tissues"]) & set(sb["tissues"]))
        common_set = frozenset(common)
        idx_a = {t: i for i, t in enumerate(sa["tissues"])}
        idx_b = {t: i for i, t in enumerate(sb["tissues"])}

        pairs = list(combinations(common, 2))
        if len(pairs) >= 2:
            coph_a = np.array([sa["cophenetic"][idx_a[t1], idx_a[t2]] for t1, t2 in pairs])
            coph_b = np.array([sb["cophenetic"][idx_b[t1], idx_b[t2]] for t1, t2 in pairs])
            r, _ = pearsonr(coph_a, coph_b)
        else:
            r = np.nan
        dist = 1 - r

        if len(common) >= 2:
            clusters_a = [sa["cluster_of"][t] for t in common]
            clusters_b = [sb["cluster_of"][t] for t in common]
            ari = adjusted_rand_score(clusters_a, clusters_b)
        else:
            ari = np.nan

        splits_a = restricted_splits(sa["splits"], common_set)
        splits_b = restricted_splits(sb["splits"], common_set)
        union = splits_a | splits_b
        rf_count = len(splits_a ^ splits_b)
        rf_jaccard = (rf_count / len(union)) if union else np.nan

        dist_matrix.loc[a, b] = dist_matrix.loc[b, a] = dist
        rf_matrix.loc[a, b] = rf_matrix.loc[b, a] = rf_jaccard
        long_rows.append({
            "scenario_a": a, "scenario_b": b, "n_common_tissues": len(common),
            "cophenetic_pearson_r": r, "distance_1_minus_r": dist, "ari": ari,
            "n_splits_a": len(splits_a), "n_splits_b": len(splits_b),
            "rf_count": rf_count, "rf_jaccard_distance": rf_jaccard,
        })

    return dist_matrix, rf_matrix, pd.DataFrame(long_rows)
