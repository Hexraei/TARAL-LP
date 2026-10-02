"""unbounded_relax: unbounded LP relaxation: integer-feasible ray (unbounded) versus no integer point (infeasible)."""

from base import Case
from mpsio import Mps

NAME = "unbounded_relax"
COUNT = 40
SEED_BASE = 56000
MIP = True


def gen(rng, seed):
    """Unbounded relaxation: integer-feasible with a ray (unbounded) or no integer point (infeasible)."""
    n = rng.randint(3, 5)
    mp = Mps("UBR")
    mp.row("N", "OBJ")
    mp.row("G", "R0")
    for j in range(n):
        c = "X%d" % j
        mp.col(c, integer=(j < n - 1))
        mp.add(c, "OBJ", -1 if j == 0 else rng.randint(0, 2))
        mp.add(c, "R0", 1 if j == 0 else 0)
    infeasible = seed % 2 == 1
    if infeasible:  # 2 x1 = 1 forces a non-integer value for x1 (integer) while x0 is a free ray
        mp.row("E", "R1")
        mp.add("X1", "R1", 2)
        mp.rhs.append(("R1", 1))
    mp.rhs.append(("R0", 0))
    return Case("unbounded_relax", seed, mp.render(), mip=True, expect="infeasible" if infeasible else "unbounded")
