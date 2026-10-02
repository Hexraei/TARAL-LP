"""general_int: general integers with box/lower bounds via MARKER blocks or LI/UI bound types."""

from base import Case, dims
from core import rand_core, encode
from mpsio import INF

NAME = "general_int"
COUNT = 80
SEED_BASE = 47000
MIP = True


def _mip_core(rng, m, n, **kw):
    kw.setdefault("integer", 0.7)
    return rand_core(rng, m, n, **kw)


def gen(rng, seed):
    m, n = dims(rng, 2, 8)
    c = _mip_core(rng, m, n, colkinds=("box", "nonneg", "lowneg"), integer=1.0)
    st = {"int_bounds": ["marker", "special"][seed % 2]}
    for j in c.ints:
        if c.col_up[j] == INF:
            c.col_up[j] = c.x0[j] + rng.randint(2, 8)
    return Case("general_int", seed, encode(c, rng, st).render(), mip=True, expect="optimal")
