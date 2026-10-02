"""big_values: MILP with all coefficients, bounds and costs scaled by 1e3..1e7."""

from base import Case, dims
from core import rand_core, encode
from mpsio import INF

NAME = "big_values"
COUNT = 40
SEED_BASE = 57000
MIP = True


def _mip_core(rng, m, n, **kw):
    kw.setdefault("integer", 0.7)
    return rand_core(rng, m, n, **kw)


def gen(rng, seed):
    m, n = dims(rng, 2, 6)
    c = _mip_core(rng, m, n, integer=0.7, colkinds=("box",))
    sc = 10.0 ** rng.randint(3, 7)
    for i in range(m):
        c.A[i] = {j: v * sc for j, v in c.A[i].items()}
        c.row_lo[i] = c.row_lo[i] * sc if c.row_lo[i] > -INF else -INF
        c.row_up[i] = c.row_up[i] * sc if c.row_up[i] < INF else INF
    c.cost = [x * sc for x in c.cost]
    return Case("big_values", seed, encode(c, rng, {}).render(), mip=True, expect="optimal")
