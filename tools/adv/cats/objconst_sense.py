"""objconst_sense: objective constants (RHS on the objective row: +-7, 1e6, 0.001, ...) and OBJSENSE MAX/MIN in section, inline and MAXIMIZE spellings."""

from base import Case, dims
from core import rand_core, encode

NAME = "objconst_sense"
COUNT = 120
SEED_BASE = 40000
MIP = False


def gen(rng, seed):
    m, n = dims(rng)
    c = rand_core(rng, m, n)
    c.const = float(rng.choice([0, 7, -7, 1e6, -1e6, 0.001, 123456.789, -3.25]))
    c.maximize = seed % 2 == 0
    if c.maximize:
        c.cost = [-x for x in c.cost]
    st = {"sense": ["section", "inline"][seed % 4 // 2], "explicit_min": not c.maximize and seed % 3 == 0}
    return Case("objconst_sense", seed, encode(c, rng, st).render(), expect="optimal")
