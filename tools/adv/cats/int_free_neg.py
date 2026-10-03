"""int_free_neg: integers with free / negative / upper-negative bounds made finite."""

from base import Case, dims
from core import rand_core, encode
from mpsio import INF

NAME = "int_free_neg"
COUNT = 60
SEED_BASE = 51000
MIP = True


def _mip_core(rng, m, n, **kw):
    kw.setdefault("integer", 0.7)
    return rand_core(rng, m, n, **kw)


def gen(rng, seed):
    m, n = dims(rng, 2, 7)
    c = _mip_core(rng, m, n, colkinds=("free", "lowneg", "upneg", "box"), integer=0.8)
    for j in c.ints:  # keep the integer part bounded so B&B terminates
        if c.col_lo[j] == -INF:
            c.col_lo[j] = c.x0[j] - rng.randint(1, 6)
        if c.col_up[j] == INF:
            c.col_up[j] = c.x0[j] + rng.randint(1, 6)
    return Case("int_free_neg", seed, encode(c, rng, {}).render(), mip=True, expect="optimal")
