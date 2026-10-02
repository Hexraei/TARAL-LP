"""setcover: set covering with binary columns, up to 45 sets x 25 elements."""

from base import Case
from mpsio import Mps

NAME = "setcover"
COUNT = 60
SEED_BASE = 48000
MIP = True


def gen(rng, seed):
    n = rng.randint(5, 22) if seed % 4 else rng.randint(30, 45)
    k = rng.randint(3, 14) if seed % 4 else rng.randint(15, 25)
    mp = Mps("COVER")
    mp.row("N", "OBJ")
    for i in range(k):
        mp.row("G", "E%d" % i)
    cov = [[i for i in range(k) if rng.random() < 0.35] for _ in range(n)]
    for i in range(k):  # every element is coverable
        if not any(i in s for s in cov):
            cov[rng.randrange(n)].append(i)
    for j in range(n):
        c = "S%d" % j
        mp.col(c, integer=True)
        mp.add(c, "OBJ", rng.randint(1, 9))
        for i in sorted(set(cov[j])):
            mp.add(c, "E%d" % i, 1)
    for i in range(k):
        mp.rhs.append(("E%d" % i, 1))
    for j in range(n):
        mp.bounds.append(("UP", "S%d" % j, 1))
    return Case("setcover", seed, mp.render(), mip=True, expect="optimal")
