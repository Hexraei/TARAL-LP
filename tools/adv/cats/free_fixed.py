"""free_fixed: free (FR) and fixed (FX or LO=UP) variables, incl. free columns with no entries and fixed columns that shift row activities."""

from base import Case, dims
from core import rand_core, encode
from mpsio import INF

NAME = "free_fixed"
COUNT = 150
SEED_BASE = 30000
MIP = False


def gen(rng, seed):
    m, n = dims(rng)
    kinds = [("free", "fixed", "nonneg", "box"), ("free", "free", "nonneg"), ("fixed", "fixed", "box", "nonneg")][seed % 3]
    c = rand_core(rng, m, n, colkinds=kinds, dens=0.5)
    st = {}
    if seed % 5 == 0:  # free column with no entries and zero cost / a fixed column with no entries
        c.n += 2
        c.col_lo += [-INF, 3.0]
        c.col_up += [INF, 3.0]
        c.x0 += [0.0, 3.0]
        c.cost += [0.0, 2.0]
    return Case("free_fixed", seed, encode(c, rng, st).render(), expect="optimal")
