"""huge_bounds: bounds, RHS and RANGES written as 1e20/1e30 (infinity by common convention), huge finite bounds (1e15), an unbounded ray bounded only by the sentinel."""

from base import Case, dims
from core import rand_core, encode
from mpsio import INF

NAME = "huge_bounds"
COUNT = 150
SEED_BASE = 31000
MIP = False


def gen(rng, seed):
    m, n = dims(rng)
    v = seed % 4
    c = rand_core(rng, m, n, colkinds=("nonneg", "free", "lowneg", "upneg", "box"))
    sp = ["1e20", "1e30", "inf"][seed % 3]
    if v == 2:  # finite but enormous bounds the optimum may sit at (not infinity)
        j = rng.randrange(n)
        if c.col_up[j] == INF:
            c.col_up[j] = 1e15
        if c.col_lo[j] == -INF:
            c.col_lo[j] = -1e15
    mp = encode(c, rng, {"huge": sp})
    if v == 1:  # rows whose RHS is the infinite sentinel: always satisfied
        mp.row("L", "RINF")
        mp.row("G", "RNINF")
        for r in ("RINF", "RNINF"):
            for j in range(min(n, 3)):
                mp.add("C%d" % j, r, 1)
        big = {"1e20": 1e20, "1e30": 1e30, "inf": 1e30}[sp]
        mp.rhs.append(("RINF", "1e+30" if sp == "inf" else big))
        mp.rhs.append(("RNINF", -big))
        return Case("huge_bounds", seed, mp.render(), expect="optimal")
    if v == 3:  # RANGES holding the sentinel on L/G rows: that side of the row is infinite (no expect: may unbound)
        have = {r for r, _ in mp.ranges}
        for i in range(c.m):
            if mp.rows[i + 1][0] in "LG" and "R%d" % i not in have and rng.random() < 0.5:
                mp.ranges.append(("R%d" % i, 1e30))
        return Case("huge_bounds", seed, mp.render())
    if v == 0 and seed % 8 == 0:  # unbounded ray whose upper bound is the sentinel: infinite, so the LP is unbounded
        mp.col("CRAY")
        mp.add("CRAY", "OBJ", -1)
        mp.bounds.append(("UP", "CRAY", {"1e20": 1e20, "1e30": 1e30, "inf": "1e+30"}[sp]))
        return Case("huge_bounds", seed, mp.render(), expect="unbounded")
    return Case("huge_bounds", seed, mp.render(), expect="optimal")
