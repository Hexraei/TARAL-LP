"""degenerate: assignment and transportation LPs with heavy cost ties, many constraints tight at one vertex (redundant combinations of tight rows) and Klee-Minty cubes: primal/dual degeneracy and stalling."""

from base import Case
from core import Core, rand_core, encode
from mpsio import Mps

NAME = "degenerate"
COUNT = 200
SEED_BASE = 26000
MIP = False


def gen(rng, seed):
    v = seed % 4
    if v == 0:  # assignment problem with heavy ties (degenerate, N=2..7)
        n = rng.randint(2, 7)
        mp = Mps("ASSIGN")
        mp.row("N", "OBJ")
        for i in range(n):
            mp.row("E", "A%d" % i)
            mp.row("E", "B%d" % i)
        for i in range(n):
            for j in range(n):
                c = "X%d_%d" % (i, j)
                mp.add(c, "OBJ", rng.choice([0, 1, 1, 2]))
                mp.add(c, "A%d" % i, 1)
                mp.add(c, "B%d" % j, 1)
        for i in range(n):
            mp.rhs.append(("A%d" % i, 1))
            mp.rhs.append(("B%d" % i, 1))
        return Case("degenerate", seed, mp.render(), expect="optimal")
    if v == 1:  # transportation with equal supplies/demands (primal degenerate)
        s, d = rng.randint(2, 5), rng.randint(2, 5)
        tot = s * d
        mp = Mps("TRANSP")
        mp.row("N", "OBJ")
        for i in range(s):
            mp.row("L", "S%d" % i)
        for j in range(d):
            mp.row("G", "D%d" % j)
        for i in range(s):
            for j in range(d):
                c = "X%d_%d" % (i, j)
                mp.add(c, "OBJ", rng.choice([1, 2, 2, 3]))
                mp.add(c, "S%d" % i, 1)
                mp.add(c, "D%d" % j, 1)
        for i in range(s):
            mp.rhs.append(("S%d" % i, d))
        for j in range(d):
            mp.rhs.append(("D%d" % j, s))
        return Case("degenerate", seed, mp.render(), expect="optimal")
    if v == 2:  # many constraints tight at one vertex: redundant combinations of tight rows
        n = rng.randint(2, 6)
        c = rand_core(rng, 1, n, dens=1.0, colkinds=("nonneg",), rowkinds=("le",), slack0=1.0)
        x0 = [float(rng.choice([0, 0, 1, 2])) for _ in range(n)]
        base = []
        for _ in range(n):
            base.append({j: float(rng.randint(-3, 3)) for j in range(n)})
        rows = list(base)
        for _ in range(rng.randint(2, 6)):  # nonneg integer combinations of tight rows stay tight at x0
            lam = [rng.randint(0, 2) for _ in base]
            rows.append({j: sum(l * r[j] for l, r in zip(lam, base)) for j in range(n)})
        cc = Core(len(rows), n)
        cc.x0 = x0
        for i, r in enumerate(rows):
            cc.A[i] = {j: val for j, val in r.items() if val != 0}
            cc.row_up[i] = sum(val * x0[j] for j, val in cc.A[i].items())
        cc.col_up = [x0[j] + rng.randint(1, 6) for j in range(n)]  # box makes it bounded
        cc.cost = [float(rng.randint(-4, 4)) for _ in range(n)]
        return Case("degenerate", seed, encode(cc, rng, {"split": 0}).render(), expect="optimal")
    # klee-minty style cube (n=3..9), scaled randomly by powers of two
    n = rng.randint(3, 9)
    mp = Mps("KM")
    mp.row("N", "OBJ")
    for i in range(n):
        mp.row("L", "R%d" % i)
    for j in range(n):
        c = "X%d" % j
        mp.add(c, "OBJ", -(2 ** (n - 1 - j)))
        for i in range(j + 1, n):
            mp.add(c, "R%d" % i, 2 ** (i - j + 1))
        mp.add(c, "R%d" % j, 1)
    for i in range(n):
        mp.rhs.append(("R%d" % i, 5 ** (i + 1)))
    return Case("degenerate", seed, mp.render(), expect="optimal")
