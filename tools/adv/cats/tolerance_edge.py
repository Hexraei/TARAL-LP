"""tolerance_edge: infeasible by exactly 1e-2..1e-5 (a row demanding the box maximum plus delta) and the delta=0 feasible twin."""

from base import Case, dims
from core import rand_core, encode
from mpsio import INF

NAME = "tolerance_edge"
COUNT = 100
SEED_BASE = 44000
MIP = False


def gen(rng, seed):
    """Boxed core plus a row demanding sum >= (max possible) + delta: infeasible by exactly delta (1e-2..1e-5);
    or the same row with delta = 0 (feasible, a single point)."""
    m, n = dims(rng, 2, 7)
    c = rand_core(rng, m, n, colkinds=("box",))
    row = {j: float(rng.choice([1, 2, 3])) for j in range(n)}
    top = sum(a * c.col_up[j] for j, a in row.items())
    delta = rng.choice([1e-2, 1e-3, 1e-4, 1e-5]) if seed % 4 else 0.0
    if not delta:  # feasible twin: the only point is x = col_up, so re-centre every core row on that point
        for i in range(c.m):
            sh = sum(v * (c.col_up[j] - c.x0[j]) for j, v in c.A[i].items())
            c.row_lo[i] += sh
            c.row_up[i] += sh
        c.x0 = list(c.col_up)
    c.A.append(row)
    c.row_lo.append(top + delta), c.row_up.append(INF)
    c.m += 1
    return Case("tolerance_edge", seed, encode(c, rng, {}).render(), expect="infeasible" if delta else "optimal")
