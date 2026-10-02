"""equality_int: integer equality systems."""

from base import Case, dims
from core import rand_core, encode

NAME = "equality_int"
COUNT = 60
SEED_BASE = 55000
MIP = True


def _mip_core(rng, m, n, **kw):
    kw.setdefault("integer", 0.7)
    return rand_core(rng, m, n, **kw)


def gen(rng, seed):
    m, n = dims(rng, 2, 6)
    c = _mip_core(rng, m, n, rowkinds=("eq", "eq", "le"), integer=1.0, colkinds=("box",), dens=0.6)
    return Case("equality_int", seed, encode(c, rng, {}).render(), mip=True, expect="optimal")
