"""2x2 histogram grid of per-gene (multi-transcript genes only, see
compute_multitranscript_variance.py) tau_range, tau_var (linear, tau in [0,1]) and
max_tpm_range, max_tpm_var (log10(x+1) -- TPM spans orders of magnitude, and x+1
handles the frequent exact-0 case, e.g. transcripts of a gene that are all equally
unexpressed)."""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IN_CSV = sys.argv[1]
OUT_PNG = sys.argv[2]

df = pd.read_csv(IN_CSV)

fig, axes = plt.subplots(2, 2, figsize=(10, 8))
panels = [
    (df["tau_range"].dropna(), "tau_range", False),
    (df["tau_var"].dropna(), "tau_var", False),
    (np.log10(df["max_tpm_range"] + 1), "log10(max_tpm_range + 1)", True),
    (np.log10(df["max_tpm_var"] + 1), "log10(max_tpm_var + 1)", True),
]
for ax, (vals, label, _) in zip(axes.flat, panels):
    ax.hist(vals, bins=30, color="steelblue", edgecolor="black")
    ax.set_xlabel(label)
    ax.set_ylabel("multi-transcript genes")

fig.suptitle(f"gtex_txs: tau and max-TPM spread among a gene's transcripts (n={len(df)} genes)")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)

print(f"genes={len(df)} -> {OUT_PNG}", file=sys.stderr)
