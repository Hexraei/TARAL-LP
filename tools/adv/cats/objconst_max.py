"""objconst_max: MILP with objective constants and OBJSENSE MAX."""

from base import Case, dims
from core import rand_core, encode
from mpsio import INF

NAME = "objconst_max"
COUNT = 50
SEED_BASE = 54000
MIP = True


def _mip_core(rng, m, n, **kw):
    kw.setdefault("integer", 0.7)
    return rand_core(rng, m, n, **kw)


def gen(rng, seed):
    m, n = dims(rng, 2, 8)
    c = _mip_core(rng, m, n, integer=0.7, colkinds=("box", "nonneg"))
    for j in c.ints:
        if c.col_up[j] == INF:
            c.col_up[j] = c.x0[j] + rng.randint(2, 8)
    c.const = float(rng.choice([5, -5, 1e5, 0.5]))
    c.maximize = True
    c.cost = [-x for x in c.cost]
    return Case("objconst_max", seed, encode(c, rng, {"sense": ["section", "inline"][seed % 2]}).render(), mip=True,
                expect="optimal")
