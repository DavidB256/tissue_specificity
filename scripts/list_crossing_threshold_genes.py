"""List genes whose tau crosses the 0.5 or 0.8 threshold going from tau_unweighted
(full resolution) to tau_free_weighted (fit to a leg's collapsed reference) -- same
categorization as plot_tau_scenarios_scatter_thresholds.py, minus the "no threshold
crossed" catch-all. Writes 4 .txt files (one per category, one gene per line, sorted
by descending |tau_free_weighted - tau_unweighted|).
Args: IN_WEIGHTED_TAU_CSV OUT_DIR PREFIX"""
import sys
import os
import pandas as pd

IN_WEIGHTED_TAU_CSV = sys.argv[1]
OUT_DIR = sys.argv[2]
PREFIX = sys.argv[3]

df = pd.read_csv(IN_WEIGHTED_TAU_CSV)[["gene", "tau_unweighted", "tau_free_weighted"]].dropna()
u, w = df["tau_unweighted"], df["tau_free_weighted"]


def categorize(u, w):
    if u > 0.8 and w < 0.8:
        return "above0.8_to_below"
    if u < 0.8 and w > 0.8:
        return "below0.8_to_above"
    if u > 0.5 and w < 0.5:
        return "above0.5_to_below"
    if u < 0.5 and w > 0.5:
        return "below0.5_to_above"
    return None


df["category"] = [categorize(ui, wi) for ui, wi in zip(u, w)]
df["abs_delta"] = (w - u).abs()

os.makedirs(OUT_DIR, exist_ok=True)
for category in ("above0.5_to_below", "below0.5_to_above", "above0.8_to_below", "below0.8_to_above"):
    sub = df[df["category"] == category].sort_values("abs_delta", ascending=False)
    out_path = f"{OUT_DIR}/{PREFIX}_{category}.txt"
    with open(out_path, "w") as f:
        f.write("\n".join(sub["gene"]) + "\n")
    print(f"{category}: n={len(sub)} -> {out_path}", file=sys.stderr)
