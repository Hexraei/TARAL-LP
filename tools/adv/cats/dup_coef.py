"""dup_coef: the same (row, column) coefficient listed two or three times in COLUMNS (sums, pairs that cancel to zero, duplicate objective entries)."""

from base import Case, dims
from core import rand_core, encode

NAME = "dup_coef"
COUNT = 150
SEED_BASE = 28000
MIP = False


def gen(rng, seed):
    m, n = dims(rng)
    c = rand_core(rng, m, n, coef="dyadic", colkinds=("nonneg", "box", "lowneg"))
    st = {"split": rng.choice([0.3, 0.6, 1.0]), "coef": "dyadic"}
    mp = encode(c, rng, st)
    if rng.random() < 0.4:  # a pair that cancels exactly: duplicate entries sum to zero
        j = rng.randrange(n)
        r = "R%d" % rng.randrange(m)
        a = float(rng.randint(1, 4))
        mp.add("C%d" % j, r, a)
        mp.add("C%d" % j, r, -a)
    if rng.random() < 0.3:  # duplicate objective coefficient
        j = rng.randrange(n)
        mp.add("C%d" % j, "OBJ", 0.5)
        mp.add("C%d" % j, "OBJ", -0.5)
    return Case("dup_coef", seed, mp.render(), expect="optimal")
