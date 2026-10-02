"""fixed_charge: fixed-charge networks with big-M links x <= M y, M up to 1e5."""

from base import Case
from mpsio import Mps

NAME = "fixed_charge"
COUNT = 60
SEED_BASE = 50000
MIP = True


def gen(rng, seed):
    n = rng.randint(2, 7)
    mp = Mps("FC")
    mp.row("N", "OBJ")
    mp.row("G", "DEM")
    M = rng.choice([50, 100, 1000, 1e5])
    for j in range(n):
        mp.row("L", "L%d" % j)
    for j in range(n):
        x, y = "X%d" % j, "Y%d" % j
        mp.add(x, "OBJ", rng.randint(1, 5))
        mp.add(x, "DEM", 1)
        mp.add(x, "L%d" % j, 1)
        mp.col(y, integer=True)
        mp.add(y, "OBJ", rng.randint(5, 60))
        mp.add(y, "L%d" % j, -M)
        mp.bounds.append(("UP", y, 1))
    mp.rhs.append(("DEM", rng.randint(5, 40)))
    return Case("fixed_charge", seed, mp.render(), mip=True, expect="optimal")
