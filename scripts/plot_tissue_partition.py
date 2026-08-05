"""Per-tissue figure: (i) the continuous exp_beta mixture fit against the whole
tissue's log10(TPM), (ii) the grid-search 'exp' group's data against its fitted
exponential, (iii) the grid-search 'beta' group's data against its fitted beta.
Subtitles report each component's share of the data, and (for ii/iii) a description
of the discriminant line that produced the hard partition."""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mixture_core import logpdf
from model_registry import MODELS

TISSUE = sys.argv[1]
IN_PARQUET = sys.argv[2]
IN_EXP_BETA_FITS = sys.argv[3]
IN_PARTITION_SUMMARY = sys.argv[4]
IN_PARTITION_LABELS = sys.argv[5]
OUT_PNG = sys.argv[6]

COMP_TYPES = MODELS["exp_beta"]


def discriminant_description(theta_deg, offset):
    theta = np.radians(theta_deg)
    w_expr, w_tau = np.cos(theta), np.sin(theta)
    return (f"boundary: {w_expr:.2f}*rank(log10TPM) + {w_tau:.2f}*rank(tau) >= {offset:.2f} "
            f"(angle={theta_deg:.0f} deg)")


long_df = pd.read_parquet(IN_PARQUET)
x_all = long_df.loc[long_df.tissue == TISSUE, ["gene", "log10_tpm"]]

cont = pd.read_csv(IN_EXP_BETA_FITS).set_index("tissue").loc[TISSUE]
part = pd.read_csv(IN_PARTITION_SUMMARY).set_index("tissue").loc[TISSUE]
labels = pd.read_parquet(IN_PARTITION_LABELS)
labels = labels.loc[labels.tissue == TISSUE, ["gene", "group"]]

merged = x_all.merge(labels, on="gene", how="inner")
x_exp_group = merged.loc[merged.group == "exp", "log10_tpm"].values
x_beta_group = merged.loc[merged.group == "beta", "log10_tpm"].values
x = x_all["log10_tpm"].values

desc = discriminant_description(part["theta_deg"], part["offset"])

fig, axes = plt.subplots(1, 3, figsize=(18, 5))

# (i) continuous exp_beta fit, whole data
ax = axes[0]
ax.hist(x, bins=60, density=True, alpha=0.5, color="steelblue")
xs = np.linspace(x.min(), x.max(), 500)
weights = [cont["weight_1"], cont["weight_2"]]
params = [{"rate": cont["comp1_rate"]},
          {"a": cont["comp2_a"], "b": cont["comp2_b"], "loc": cont["comp2_loc"], "scale": cont["comp2_scale"]}]
total = np.zeros_like(xs)
for ct, w, p in zip(COMP_TYPES, weights, params):
    comp = w * np.exp(logpdf(xs, ct, p))
    total += comp
    ax.plot(xs, comp, "--", lw=1.5)
ax.plot(xs, total, color="darkred", lw=2)
ax.set_title(f"Continuous exp_beta fit (whole data, n={len(x)})\n"
             f"weight: exp={weights[0]:.1%}, beta={weights[1]:.1%}", fontsize=10)
ax.set_xlabel("log10(TPM)")

# (ii) grid-search exp group only
ax = axes[1]
ax.hist(x_exp_group, bins=40, density=True, alpha=0.5, color="steelblue")
xs2 = np.linspace(x_exp_group.min(), x_exp_group.max(), 300)
ax.plot(xs2, np.exp(logpdf(xs2, "expon", {"rate": part["exp_rate"]})), color="darkred", lw=2)
frac_exp = part["n_exp"] / part["n"]
ax.set_title(f"Grid-search EXP group (n={part['n_exp']}, {frac_exp:.1%} of data)\n{desc}", fontsize=8)
ax.set_xlabel("log10(TPM)")

# (iii) grid-search beta group only
ax = axes[2]
ax.hist(x_beta_group, bins=40, density=True, alpha=0.5, color="steelblue")
xs3 = np.linspace(x_beta_group.min(), x_beta_group.max(), 300)
ax.plot(xs3, np.exp(logpdf(xs3, "beta", {"a": part["beta_a"], "b": part["beta_b"],
                                          "loc": part["beta_loc"], "scale": part["beta_scale"]})),
        color="darkred", lw=2)
frac_beta = part["n_beta"] / part["n"]
ax.set_title(f"Grid-search BETA group (n={part['n_beta']}, {frac_beta:.1%} of data)\n{desc}", fontsize=8)
ax.set_xlabel("log10(TPM)")

fig.suptitle(f"{TISSUE}: exp_beta continuous fit vs grid-search hard partition (BIC delta={part['bic_delta']:.0f})")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)
