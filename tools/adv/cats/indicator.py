"""indicator: big-M indicator links with M up to 1e7 and binary y: integrality tolerance against big-M."""

from base import Case
from mpsio import Mps

NAME = "indicator"
COUNT = 60
SEED_BASE = 58000
MIP = True


def gen(rng, seed):
    """Big-M links x <= M*y with M up to 1e7 and a demand row: the LP relaxation lets y be ~1e-7, so integrality
    tolerance handling decides correctness (y = 1e-7 'binary' would pay almost nothing for 1e7 of capacity)."""
    n = rng.randint(2, 6)
    M = rng.choice([1e4, 1e5, 1e6, 1e7])
    mp = Mps("BIGM")
    mp.row("N", "OBJ")
    mp.row("G", "DEM")
    for j in range(n):
        mp.row("L", "L%d" % j)
    cap = [rng.randint(1, 5) * 1000 for _ in range(n)]
    for j in range(n):
        x, y = "X%d" % j, "Y%d" % j
        mp.add(x, "OBJ", rng.randint(1, 4))
        mp.add(x, "DEM", 1)
        mp.add(x, "L%d" % j, 1)
        mp.col(y, integer=True)
        mp.add(y, "OBJ", rng.randint(10, 500))
        mp.add(y, "L%d" % j, -M)
        mp.bounds.append(("UP", x, cap[j]))
        mp.bounds.append(("UP", y, 1))
    mp.rhs.append(("DEM", rng.randint(500, sum(cap) // 2)))
    return Case("indicator", seed, mp.render(), mip=True, expect="optimal")
