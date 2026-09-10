"""Vertical bar chart of the free_weights per-tissue weight (winner of N DE chains,
see optimize_free_weights.py / combine_free_weights_chains.py) for every GTEx tissue,
sorted descending. Brain_* tissues highlighted separately from all other tissues.
Weight axis is strictly bounded to [0, 1]."""
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IN_WEIGHTS = sys.argv[1]
OUT_PNG = sys.argv[2]

df = pd.read_csv(IN_WEIGHTS).sort_values("weight", ascending=False)


def colors_of(d):
    return d["tissue"].str.contains("Brain").map({True: "tab:red", False: "tab:blue"})


legend_handles = [
    plt.Rectangle((0, 0), 1, 1, color="tab:red", label="Brain_*"),
    plt.Rectangle((0, 0), 1, 1, color="tab:blue", label="other"),
]
title = f"free_weights: optimized per-tissue weight (n={len(df)} tissues, winner chain)"

fig, ax = plt.subplots(figsize=(16, 8))
ax.bar(df["tissue"], df["weight"], color=colors_of(df))
ax.set_ylim(0, 1)
ax.set_ylabel("free weight")
ax.set_xlabel("tissue")
ax.tick_params(axis="x", labelsize=6, rotation=90)
ax.legend(handles=legend_handles, loc="upper right", fontsize=8)
ax.set_title(title)
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)

print(f"tissues={len(df)}", file=sys.stderr)
