"""big_status: 500-2500 row sparse LPs whose status is made infeasible or unbounded by construction (alternating):
infeasible = a feasible big core plus one row that demands more than the box allows (Farkas row on 3-6 columns, far
from the sparsity pattern the rest of the model has); unbounded = a feasible big core plus one free column with a
negative cost that only appears in a row it can drift along (a recession ray). Detecting these at scale is what the
interior-point and dual paths have to do without a basis to read a certificate from."""

from base import Case
from core import rand_core, encode
from mpsio import INF

NAME = "big_status"
COUNT = 24
SEED_BASE = 63000
MIP = False
BIG = True


def gen(rng, seed):
    m = int(round(500 * (5.0 ** rng.random())))
    n = int(m * rng.uniform(1.2, 1.8))
    c = rand_core(rng, m, n, nnz=rng.randint(4, 7), colkinds=("nonneg", "box", "lowneg"),
                  rowkinds=("eq", "le", "ge", "rng"))
    if seed % 2 == 0:  # infeasible: sum of 3..6 boxed columns must exceed its maximum
        cols = rng.sample(range(n), rng.randint(3, 6))
        for j in cols:
            c.col_lo[j], c.col_up[j] = 0.0, float(rng.randint(1, 4))
            c.x0[j] = 0.0
        mx = sum(c.col_up[j] for j in cols)
        c.A.append({j: 1.0 for j in cols})
        c.row_lo.append(mx + float(rng.randint(1, 3)))
        c.row_up.append(INF)
        c.m += 1
        expect = "infeasible"
    else:  # unbounded: a new column free to move, cost -1, in one new equality row with a free partner column
        a, b = n, n + 1
        c.n += 2
        c.col_lo += [0.0, -INF]
        c.col_up += [INF, INF]
        c.cost += [-1.0, 0.0]
        c.x0 += [0.0, 0.0]
        c.ints = set(c.ints)
        c.A.append({a: 1.0, b: -1.0})  # a - b = 0, b free, a >= 0: ray (1, 1)
        c.row_lo.append(0.0)
        c.row_up.append(0.0)
        c.m += 1
        expect = "unbounded"
    return Case("big_status", seed, encode(c, rng, {"split": 0.0}).render(), expect=expect, note="%d x %d" % (c.m, c.n))
