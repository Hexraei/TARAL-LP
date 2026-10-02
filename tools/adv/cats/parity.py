"""parity: even-coefficient equalities with odd RHS: LP-feasible but integer-infeasible (and feasible twins)."""

from base import Case
from mpsio import Mps

NAME = "parity"
COUNT = 40
SEED_BASE = 49000
MIP = True


def gen(rng, seed):
    """Equalities over even coefficients with an odd RHS: LP-feasible, integer-infeasible (or feasible if even)."""
    n = rng.randint(2, 8)
    mp = Mps("PAR")
    mp.row("N", "OBJ")
    k = rng.randint(1, 2)
    for r in range(k):
        mp.row("E", "R%d" % r)
    bad = seed % 3 != 0
    for j in range(n):
        c = "X%d" % j
        mp.col(c, integer=True)
        mp.add(c, "OBJ", rng.randint(-3, 3))
        for r in range(k):
            mp.add(c, "R%d" % r, 2 * rng.randint(1, 4))
    for r in range(k):
        mp.rhs.append(("R%d" % r, 2 * rng.randint(2, 20) + (1 if bad and r == 0 else 0)))
    for j in range(n):
        mp.bounds.append(("UP", "X%d" % j, 20))
    return Case("parity", seed, mp.render(), mip=True, expect="infeasible" if bad else None)
