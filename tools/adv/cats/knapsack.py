"""knapsack: binary knapsacks (BV bounds) with 1-3 rows, up to 40 items."""

from base import Case
from mpsio import Mps

NAME = "knapsack"
COUNT = 80
SEED_BASE = 46000
MIP = True


def gen(rng, seed):
    n = rng.randint(4, 20) if seed % 4 else rng.randint(25, 40)
    mp = Mps("KNAP")
    mp.row("N", "OBJ")
    k = rng.randint(1, 3)
    w = [[rng.randint(1, 30) for _ in range(n)] for _ in range(k)]
    for r in range(k):
        mp.row("L", "W%d" % r)
    v = [rng.randint(1, 40) for _ in range(n)]
    for j in range(n):
        c = "X%d" % j
        mp.col(c, integer=True)
        mp.add(c, "OBJ", -v[j])
        for r in range(k):
            mp.add(c, "W%d" % r, w[r][j])
    for r in range(k):
        mp.rhs.append(("W%d" % r, sum(w[r]) // rng.choice([2, 3, 4])))
    for j in range(n):
        mp.bounds.append(("BV", "X%d" % j, None))
    return Case("knapsack", seed, mp.render(), mip=True, expect="optimal")
