"""int_bound_types: BV, LI, UI, MI+UI and fractional LI/UI values (ceil/floor semantics) on integer columns."""

from base import Case
from core import rand_core, encode

NAME = "int_bound_types"
COUNT = 80
SEED_BASE = 53000
MIP = True


def _mip_core(rng, m, n, **kw):
    kw.setdefault("integer", 0.7)
    return rand_core(rng, m, n, **kw)


def gen(rng, seed):
    """BV / LI / UI / MI+UI on integer columns, fractional bounds that round (LI 0.5 -> 1), FX on integers."""
    n = rng.randint(3, 9)
    m = rng.randint(2, 5)
    c = _mip_core(rng, m, n, integer=1.0, colkinds=("box",))
    mp = encode(c, rng, {"int_bounds": "special"})
    j = rng.randrange(n)  # fractional LI/UI: ceil / floor semantics must change the optimum consistently
    for t in range(len(mp.bounds)):
        if mp.bounds[t][1] == "C%d" % j and mp.bounds[t][0] == "UI":
            mp.bounds[t] = ("UI", "C%d" % j, mp.bounds[t][2] + 0.5)
        if mp.bounds[t][1] == "C%d" % ((j + 1) % n) and mp.bounds[t][0] == "LI":
            mp.bounds[t] = ("LI", mp.bounds[t][1], mp.bounds[t][2] - 0.5)
    return Case("int_bound_types", seed, mp.render(), mip=True, expect="optimal")
