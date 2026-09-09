"""Bar chart of the free_weights per-tissue weight (winner of N DE chains, see
optimize_free_weights.py / combine_free_weights_chains.py) for every GTEx tissue.
Brain_* tissues highlighted separately from all other tissues. Weight axis is
strictly bounded to [0, 1]. Writes both a horizontal (tissues on y-axis, sorted
ascending) and a vertical (tissues on x-axis, sorted descending) version."""
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IN_WEIGHTS = sys.argv[1]
OUT_PNG_H = sys.argv[2]
OUT_PNG_V = sys.argv[3]

df = pd.read_csv(IN_WEIGHTS)


def colors_of(d):
    return d["tissue"].str.contains("Brain").map({True: "tab:red", False: "tab:blue"})


legend_handles = [
    plt.Rectangle((0, 0), 1, 1, color="tab:red", label="Brain_*"),
    plt.Rectangle((0, 0), 1, 1, color="tab:blue", label="other"),
]
title = f"free_weights: optimized per-tissue weight (n={len(df)} tissues, winner chain)"

df_h = df.sort_values("weight", ascending=True)
fig, ax = plt.subplots(figsize=(8, 16))
ax.barh(df_h["tissue"], df_h["weight"], color=colors_of(df_h))
ax.set_xlim(0, 1)
ax.set_xlabel("free weight")
ax.set_ylabel("tissue")
ax.tick_params(axis="y", labelsize=6)
ax.legend(handles=legend_handles, loc="lower right", fontsize=8)
ax.set_title(title)
fig.tight_layout()
fig.savefig(OUT_PNG_H, dpi=300)
plt.close(fig)

df_v = df.sort_values("weight", ascending=False)
fig, ax = plt.subplots(figsize=(16, 8))
ax.bar(df_v["tissue"], df_v["weight"], color=colors_of(df_v))
ax.set_ylim(0, 1)
ax.set_ylabel("free weight")
ax.set_xlabel("tissue")
ax.tick_params(axis="x", labelsize=6, rotation=90)
ax.legend(handles=legend_handles, loc="upper right", fontsize=8)
ax.set_title(title)
fig.tight_layout()
fig.savefig(OUT_PNG_V, dpi=300)
plt.close(fig)

print(f"tissues={len(df)}", file=sys.stderr)
