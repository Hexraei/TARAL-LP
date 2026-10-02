"""big_assign: wide degenerate LPs: assignment (N = 25..70, so up to 140 rows x 4900 columns) and transportation problems
with heavy cost ties and equal supplies/demands. Many more columns than rows, every basis heavily degenerate."""

from base import Case
from mpsio import Mps

NAME = "big_assign"
COUNT = 16
SEED_BASE = 62000
MIP = False
BIG = True


def gen(rng, seed):
    mp = Mps("WIDE")
    mp.row("N", "OBJ")
    if seed % 2 == 0:
        n = rng.randint(25, 70)
        for i in range(n):
            mp.row("E", "A%d" % i)
            mp.row("E", "B%d" % i)
        for i in range(n):
            for j in range(n):
                c = "X%d_%d" % (i, j)
                mp.add(c, "OBJ", rng.choice([0, 1, 1, 2, 3]))
                mp.add(c, "A%d" % i, 1)
                mp.add(c, "B%d" % j, 1)
        for i in range(n):
            mp.rhs += [("A%d" % i, 1), ("B%d" % i, 1)]
        note = "assignment %d" % n
    else:
        s, d = rng.randint(20, 60), rng.randint(20, 80)
        sup = [rng.randint(5, 9) for _ in range(s)]
        dem = [0] * d
        for _ in range(sum(sup)):  # demands sum exactly to supplies
            dem[rng.randrange(d)] += 1
        for i in range(s):
            mp.row("L", "S%d" % i)
        for j in range(d):
            mp.row("G", "D%d" % j)
        for i in range(s):
            for j in range(d):
                c = "X%d_%d" % (i, j)
                mp.add(c, "OBJ", rng.choice([1, 2, 2, 3, 4]))
                mp.add(c, "S%d" % i, 1)
                mp.add(c, "D%d" % j, 1)
        for i in range(s):
            mp.rhs.append(("S%d" % i, sup[i]))
        for j in range(d):
            mp.rhs.append(("D%d" % j, dem[j]))
        note = "transportation %dx%d" % (s, d)
    return Case("big_assign", seed, mp.render(), expect="optimal", note=note)
