"""negzero: '-0' spelled in coefficients, RHS, RANGES, bounds, the objective constant and costs, incl. UP -0 on a column with the default lower bound."""

from base import Case, dims
from core import rand_core, encode

NAME = "negzero"
COUNT = 60
SEED_BASE = 32000
MIP = False


def gen(rng, seed):
    m, n = dims(rng)
    c = rand_core(rng, m, n, colkinds=("nonneg", "box", "free", "fixed", "upneg"))
    if seed % 2:
        c.const = 0.0
    mp = encode(c, rng, {"negzero": True, "zeros": 3, "zero_const_line": True, "zero_rhs_lines": True,
                         "zero_cost_lines": True, "explicit_min": seed % 3 == 0})
    mp.bounds = [(t, col, "-0" if (v is not None and not isinstance(v, str) and v == 0) else v) for t, col, v in mp.bounds]
    if seed % 4 == 0:  # UP -0 on a column with the default lower bound 0: ub = 0, lb stays 0 (not -inf)
        mp.col("NZ")
        mp.add("NZ", "OBJ", -1)
        mp.bounds.append(("UP", "NZ", "-0"))
    return Case("negzero", seed, mp.render(), expect="optimal")
