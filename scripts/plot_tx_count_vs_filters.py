"""Number of transcripts retained vs. two alternative expression filters, on one graph
with two x-axes (bottom: median TPM threshold; top: donor-fraction threshold). A
transcript is retained if, in at least one tissue, (a) median TPM >= T (T=1..20; GCT is
already floored at TPM<1), or (b) TPM > 0 in >= X% of that tissue's donors (X=0..100,
step 10; X=0 keeps everything)."""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IN_GCT, IN_DONOR_FRAC_CSV, OUT_PNG, OUT_CSV = sys.argv[1:5]

gct = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
med_max = gct[list(gct.columns[2:])].max(axis=1).values
frac = pd.read_csv(IN_DONOR_FRAC_CSV).set_index("transcript_id")
frac_max = frac.max(axis=1).reindex(gct["Name"]).values

tpm_thr = np.arange(1, 21)
don_thr = np.arange(0, 101, 10)
n_tpm = np.array([(med_max >= t).sum() for t in tpm_thr])
n_don = np.array([(frac_max >= x / 100 - 1e-12).sum() for x in don_thr])
pd.concat([pd.DataFrame({"filter": "median_tpm", "threshold": tpm_thr, "n_tx": n_tpm}),
           pd.DataFrame({"filter": "donor_pct", "threshold": don_thr, "n_tx": n_don})]).to_csv(OUT_CSV, index=False)

fig, ax = plt.subplots(figsize=(8, 5))
c1, c2 = "tab:blue", "tab:orange"
ax.plot(tpm_thr, n_tpm, marker="o", color=c1)
ax.set_xlabel("median TPM threshold (>= T in at least one tissue)", color=c1)
ax.set_xticks(tpm_thr)
ax.tick_params(axis="x", colors=c1)
ax.set_ylabel("number of transcripts")
ax.set_ylim(bottom=0)
ax2 = ax.twiny()
ax2.plot(don_thr, n_don, marker="s", color=c2)
ax2.set_xlabel("donor threshold (TPM > 0 in >= X% of donors, in at least one tissue)", color=c2)
ax2.set_xticks(don_thr)
ax2.tick_params(axis="x", colors=c2)
ax.set_title("gtex_v11_tx_protein_coding: transcripts retained", y=1.15)
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)
