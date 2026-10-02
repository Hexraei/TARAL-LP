"""markers_mixed: mixed integer/continuous columns, interleaved MARKER blocks, duplicate coefficients, BV/LI/UI versus markers."""

from base import Case, dims
from core import rand_core, encode
from mpsio import INF

NAME = "markers_mixed"
COUNT = 80
SEED_BASE = 52000
MIP = True


def _mip_core(rng, m, n, **kw):
    kw.setdefault("integer", 0.7)
    return rand_core(rng, m, n, **kw)


def gen(rng, seed):
    m, n = dims(rng, 3, 10)
    c = _mip_core(rng, m, n, integer=0.5, colkinds=("nonneg", "box"))
    for j in c.ints:
        if c.col_up[j] == INF:
            c.col_up[j] = c.x0[j] + rng.randint(2, 9)
    if not c.ints:
        c.ints = {0}
        c.col_up[0] = max(c.col_up[0], c.x0[0] + 3) if c.col_up[0] < INF else c.x0[0] + 3
    st = {"split": 0.3, "int_bounds": ["marker", "special"][seed % 2]}
    return Case("markers_mixed", seed, encode(c, rng, st).render(), mip=True, expect="optimal")
