"""big_sparse: sparse LPs with 500-3000 rows (1.2-2x as many columns), 4-8 nonzeros per row, every row and column kind.
Feasible and bounded by construction (dual-built costs), so the expected status is optimal."""

from base import Case
from core import rand_core, encode

NAME = "big_sparse"
COUNT = 24
SEED_BASE = 59000
MIP = False
BIG = True


def gen(rng, seed):
    m = int(round(500 * (6.0 ** rng.random())))  # log-uniform 500..3000
    n = int(m * rng.uniform(1.2, 2.0))
    c = rand_core(rng, m, n, nnz=rng.randint(4, 8), colkinds=("nonneg", "box", "lowneg", "free"),
                  rowkinds=("eq", "le", "ge", "rng"))
    return Case("big_sparse", seed, encode(c, rng, {"split": 0.02}).render(), expect="optimal",
                note="%d x %d" % (m, n))
