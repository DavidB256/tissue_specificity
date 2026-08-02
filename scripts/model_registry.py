"""Named mixture models: abbreviation -> ordered list of component types."""
MODELS = {
    "gauss2": ["norm", "norm"],
    "gauss3": ["norm", "norm", "norm"],
    "exp_gauss": ["expon", "norm"],
    "exp_gauss2": ["expon", "norm", "norm"],
    "exp_skewgauss": ["expon", "skewnorm"],
    "exp_t": ["expon", "t"],
    "exp_gauss_t": ["expon", "norm", "t"],
    "exp_skewgauss_t": ["expon", "skewnorm", "t"],
}
