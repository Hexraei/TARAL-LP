"""ranges: RANGES on E/L/G rows (positive, negative and zero R) with mostly negative RHS, plus an N row carrying RHS and RANGES."""

from base import Case, dims
from core import rand_core, encode

NAME = "ranges"
COUNT = 200
SEED_BASE = 39000
MIP = False


def gen(rng, seed):
    """RANGES on E/L/G rows (positive, negative, zero R) with mostly negative activities/RHS, plus an N row that
    carries an RHS and a RANGES entry which must be ignored."""
    m, n = dims(rng)
    c = rand_core(rng, m, n, rowkinds=("rng", "rng", "eq", "le", "ge"), slack0=0.3,
                  colkinds=("lowneg", "upneg", "box", "free", "nonneg"))
    mp = encode(c, rng, {})
    mp.row("N", "FREE")
    for j in range(min(n, 2)):
        mp.add("C%d" % j, "FREE", rng.randint(1, 4))
    mp.rhs.append(("FREE", rng.randint(-5, 5)))
    mp.ranges.append(("FREE", rng.randint(1, 5)))
    return Case("ranges", seed, mp.render(), expect="optimal")
