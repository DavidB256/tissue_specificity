$$
\tau(\vec x) = \frac{\sum_{i=1}^{n} \left(1 - \hat{x}_{i}\right)}{n - 1}
$$

$$
\tau_{w}(\vec x;\, \vec w) = \frac{\sum_{i=1}^{n} w_i \left(1 - \hat{x}_{i}\right)}{\sum_{i=1}^{n} w_i - w_{i^*}}
$$

$$
\vec w \in (0, 1]^{N_\text{tissues}}
$$

$$
w = \Sigma^{-1}\mathbf{1}
$$

$$
\argmax_{\vec w\, \in\, (0, 1]^{N_\text{tissues}}} \rho\big(\tau_w(\vec x_\text{full};\ \vec w),\ \tau(\vec x_\text{reduced})\big)
$$