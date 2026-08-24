"""Horizontal bar chart of the free_weights per-tissue weight (winner of N DE chains,
see optimize_free_weights.py / combine_free_weights_chains.py) for every GTEx tissue,
sorted ascending. Brain_* tissues highlighted separately."""
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IN_WEIGHTS = sys.argv[1]
OUT_PNG = sys.argv[2]

df = pd.read_csv(IN_WEIGHTS).sort_values("weight", ascending=True)
is_brain = df["tissue"].str.contains("Brain")
colors = is_brain.map({True: "tab:red", False: "tab:blue"})

fig, ax = plt.subplots(figsize=(8, 16))
ax.barh(df["tissue"], df["weight"], color=colors)
ax.axvline(1 / len(df), ls="--", color="gray", lw=1, label=f"uniform (1/{len(df)})")
ax.set_xlabel("free weight")
ax.set_ylabel("tissue")
ax.tick_params(axis="y", labelsize=6)
ax.legend(handles=[
    plt.Rectangle((0, 0), 1, 1, color="tab:red", label="Brain_*"),
    plt.Rectangle((0, 0), 1, 1, color="tab:blue", label="other"),
    ax.lines[0],
], loc="lower right", fontsize=8)
ax.set_title(f"free_weights: optimized per-tissue weight (n={len(df)} tissues, winner chain)")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)
print(f"tissues={len(df)}", file=sys.stderr)
