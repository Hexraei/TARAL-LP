"""big_cover: MILP set covering with 300-1200 elements (rows) and 1.5-3x as many binary sets. Both solvers are expected to
hit their time caps on many of these; a case where the engine returns a point is checked exactly, the rest is reported
as time-outs (noans) or unresolved reference, never as a pass."""

from base import Case
from mpsio import Mps

NAME = "big_cover"
COUNT = 12
SEED_BASE = 64000
MIP = True
BIG = True


def gen(rng, seed):
    m = int(round(300 * (4.0 ** rng.random())))
    n = int(m * rng.uniform(1.5, 3.0))
    mp = Mps("COVER")
    mp.row("N", "OBJ")
    for i in range(m):
        mp.row("G", "E%d" % i)
    cover = [[] for _ in range(m)]
    for j in range(n):
        for i in rng.sample(range(m), rng.randint(2, 6)):
            cover[i].append(j)
    for i in range(m):
        if not cover[i]:
            cover[i].append(rng.randrange(n))
    for j in range(n):
        mp.col("S%d" % j, integer=True)
        mp.add("S%d" % j, "OBJ", rng.randint(1, 9))
    for i in range(m):
        for j in cover[i]:
            mp.add("S%d" % j, "E%d" % i, 1)
        mp.rhs.append(("E%d" % i, 1))
    for j in range(n):
        mp.bounds.append(("BV", "S%d" % j, None))
    return Case("big_cover", seed, mp.render(), mip=True, expect="optimal", note="%d x %d" % (m, n))
