"""large_sparse: sparse LPs with 60-150 rows and 80-200 columns mixing every row and column kind."""

from base import Case
from core import rand_core, encode

NAME = "large_sparse"
COUNT = 60
SEED_BASE = 45000
MIP = False


def gen(rng, seed):
    m, n = rng.randint(60, 150), rng.randint(80, 200)
    c = rand_core(rng, m, n, dens=rng.choice([0.02, 0.05]), colkinds=("nonneg", "box", "lowneg", "free"),
                  rowkinds=("eq", "le", "ge", "rng"))
    return Case("large_sparse", seed, encode(c, rng, {"split": 0.1}).render(), expect="optimal")
