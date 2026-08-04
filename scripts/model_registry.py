"""Named mixture models: abbreviation -> ordered list of component types."""
MODELS = {
    "gauss2": ["norm", "norm"],
    "exp_gauss": ["expon", "norm"],
    "exp_gauss2": ["expon", "norm", "norm"],
    "exp_skewgauss": ["expon", "skewnorm"],
    "exp_skewgauss2": ["expon", "skewnorm", "skewnorm"],
    "exp_skewgauss_gauss": ["expon", "skewnorm", "norm"],
    "exp_t": ["expon", "t"],
    "exp_gauss_t": ["expon", "norm", "t"],
    "exp_skewgauss_t": ["expon", "skewnorm", "t"],
    "exp_beta": ["expon", "beta"],
    "exp_gamma": ["expon", "gamma"],
}
