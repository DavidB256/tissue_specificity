# tau, tau_weighted, and the GLS weighting approach

## Tissue specificity index (tau)

For gene $g$ with expression $x_{g,1}, \dots, x_{g,n}$ across $n$ tissues (TPM,
floored: values $<1$ set to $0$), let $\hat{x}_{g,i} = x_{g,i} / \max_j x_{g,j}$.
Unweighted tau (Yanai et al. 2005):

$$
\tau_{uw} = \frac{\sum_{i=1}^{n} \left(1 - \hat{x}_{i}\right)}{n - 1}
$$

$$
\tau_{w} = \frac{\sum_{i=1}^{n} w_i \left(1 - \hat{x}_{i}\right)}{\sum_{i=1}^{n} w_i - w_{i^*}}
$$

$$
\vec w \in (0, 1]^{N_\text{tissues}}
$$

$$
w = \Sigma^{-1}\mathbf{1}
$$

$$
\argmax_{\vec w\, \in\, (0, 1]^{N_\text{clusters}}} \rho\big(\tau_w(\text{full}),\ \tau_{uw}(\text{reduced})\big)
$$

$$
\argmax_{\vec w\, \in\, (0, 1]^{N_\text{tissues}}} \rho\big(\vec\tau_w(\text{full}),\ \vec\tau_{uw}(\text{reduced})\big)
$$


$\tau_g \in [0, 1]$: $0$ = uniform expression across all tissues, $1$ = expression
confined to a single tissue.

## Weighted tau

Introduce a per-tissue weight vector $w = (w_1, \dots, w_n)$, $w_i \ge 0$. Weighted
tau replaces the uniform $1/(n-1)$ averaging with a weighted average that excludes
the max-expressing tissue $i^*_g = \arg\max_j x_{g,j}$ from the denominator:

$$
\tau_{w,g} = \frac{\sum_{i=1}^{n} w_i \left(1 - \hat{x}_{g,i}\right)}{\sum_{i=1}^{n} w_i - w_{i^*_g}}
$$

Setting $w_i = 1 \ \forall i$ recovers $\tau_g$ exactly, since the numerator's
$i=i^*_g$ term vanishes ($1-\hat x_{g,i^*_g}=0$) and the denominator becomes $n-1$.

## GLS tissue weighting

Motivation: tissues are not independent — e.g. GTEx's many `Brain_*` subregions are
highly correlated, so naive unweighted tau lets redundant tissues each pull full
weight in the specificity average, inflating apparent specificity for genes that are
merely brain-enriched. GLS weighting down-weights tissues in proportion to how
redundant they are with the rest of the panel.

Let $R \in \mathbb{R}^{n \times n}$ be the tissue $\times$ tissue Pearson correlation
matrix of $\log_{10}(\text{TPM}+1)$ across genes (restricted to genes expressed,
TPM $\ge 1$, in at least one tissue):

$$
R_{ij} = \mathrm{corr}\left(\log_{10}(x_{\cdot,i}+1),\ \log_{10}(x_{\cdot,j}+1)\right)
$$

GLS weights solve $R w = \mathbf{1}$, i.e. $w = R^{-1}\mathbf{1}$ — the standard
generalized-least-squares correction for correlated observations, applied here to
tissues rather than samples. When $R = I$ (all tissues uncorrelated), $w = \mathbf{1}$
and weighted tau reduces to unweighted tau.

Unregularized $w$ can contain negative entries under multicollinearity (reported,
not clipped, as $\tau_{w,\text{gls}}$; left-censored at $-1$ to drop unbounded
outliers). To guarantee non-negative, interpretable weights, ridge-shrink $R$:

$$
\Sigma \vec w = \mathbf1 \implies \vec w = \Sigma^{-1} \mathbf1
$$

$$
\Sigma_\lambda = (1-\lambda) \Sigma + \lambda I, \qquad w_\lambda = R_\lambda^{-1}\mathbf{1}
$$

and take the smallest $\lambda^\* \in [0,1]$ (found by binary search) such that
every entry of $w_{\lambda^\*}$ is positive:

$$
\lambda^\* = \min\left\{\lambda \in [0,1] : w_\lambda \ge 0 \text{ elementwise}\right\}
$$

$\tau_{w,\text{gls\_ridge}}$ uses $w = w_{\lambda^\*}$ in the weighted-tau formula
above. $\lambda=1$ gives $R_1 = I$, which is always all-positive ($w=\mathbf 1$), so
the search is always well-posed.
