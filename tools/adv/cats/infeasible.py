"""infeasible: structured infeasible LPs: contradictory bounds, Farkas rows, equalities without nonnegative solution, a row beyond the bound box, infeasible plus an improving ray, disjoint ranges on identical rows."""

from base import Case, dims
from core import rand_core, encode
from mpsio import Mps, INF

NAME = "infeasible"
COUNT = 150
SEED_BASE = 35000
MIP = False


def gen(rng, seed):
    v = seed % 6
    m, n = dims(rng, 2, 7)
    if v == 0:  # contradictory bounds
        c = rand_core(rng, m, n)
        j = rng.randrange(n)
        c.col_lo[j], c.col_up[j] = 5.0, 3.0
        mp = encode(c, rng, {})
        mp.bounds = [b for b in mp.bounds if b[1] != "C%d" % j]
        mp.bounds.append(("LO", "C%d" % j, 5))
        mp.bounds.append(("UP", "C%d" % j, 3))
    elif v == 1:  # sum x >= 10 and sum x <= 5
        mp = Mps("INF")
        mp.row("N", "OBJ"), mp.row("G", "A"), mp.row("L", "B")
        k = rng.randint(1, 5)
        for j in range(k):
            mp.add("X%d" % j, "OBJ", rng.randint(-2, 2))
            mp.add("X%d" % j, "A", 1)
            mp.add("X%d" % j, "B", 1)
        mp.rhs += [("A", 10), ("B", 5)]
    elif v == 2:  # equalities with no nonnegative solution
        mp = Mps("INF")
        mp.row("N", "OBJ"), mp.row("E", "A")
        k = rng.randint(1, 5)
        for j in range(k):
            mp.add("X%d" % j, "OBJ", 1)
            mp.add("X%d" % j, "A", rng.randint(1, 3))
        mp.rhs.append(("A", -rng.randint(1, 5)))
    elif v == 3:  # boxed core + a row demanding more than the box allows
        c = rand_core(rng, m, n, colkinds=("box",))
        row = {j: float(rng.randint(1, 3)) for j in range(n)}
        top = sum(a * c.col_up[j] for j, a in row.items())
        c.A.append(row)
        c.row_lo.append(top + rng.randint(1, 3)), c.row_up.append(INF)
        c.m += 1
        mp = encode(c, rng, {})
    elif v == 4:  # infeasible AND a free improving ray: primal infeasibility must win
        c = rand_core(rng, m, n, colkinds=("box",))
        row = {j: 1.0 for j in range(n)}
        c.A.append(row)
        c.row_lo.append(sum(c.col_up) + 1.0), c.row_up.append(INF)
        c.m += 1
        c.n += 1
        c.col_lo.append(0.0), c.col_up.append(INF), c.x0.append(0.0), c.cost.append(-1.0)
        mp = encode(c, rng, {})
    else:  # infeasible via a range: two rows with the same coefficients and disjoint ranges
        c = rand_core(rng, m, n)
        i = rng.randrange(m)
        c.A.append(dict(c.A[i]))
        c.m += 1
        c.row_lo.append(0.0), c.row_up.append(0.0)
        c.A.append(dict(c.A[i]))
        c.m += 1
        c.row_lo.append(1.0), c.row_up.append(2.0)
        mp = encode(c, rng, {})
    return Case("infeasible", seed, mp.render(), expect="infeasible")
