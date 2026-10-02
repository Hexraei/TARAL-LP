"""unbounded: unbounded LPs with a recession ray (row-compatible signs, free column), MIN and MAX."""

from base import Case, dims
from core import rand_core, encode
from mpsio import INF

NAME = "unbounded"
COUNT = 150
SEED_BASE = 36000
MIP = False


def gen(rng, seed):
    m, n = dims(rng, 2, 7)
    v = seed % 3
    c = rand_core(rng, m, n, colkinds=("box", "box", "nonneg") if v == 0 else ("box",))
    # one extra column k with a recession ray: row signs compatible, cost < 0
    k = n
    c.n += 1
    sign = rng.choice([1.0, -1.0])
    lo = INF
    c.col_lo.append(0.0 if v != 2 else -INF), c.col_up.append(INF), c.x0.append(0.0)
    c.cost.append(-float(rng.randint(1, 3)))
    for i in range(m):
        if rng.random() < 0.5:
            if c.row_lo[i] == -INF and c.row_up[i] < INF:
                c.A[i][k] = -float(rng.randint(0, 3)) or -1.0  # L row: increasing x_k keeps it satisfied
            elif c.row_up[i] == INF and c.row_lo[i] > -INF:
                c.A[i][k] = float(rng.randint(1, 3))
    if v == 2:  # free column: ray along -x_k needs a symmetric recession; make column entries zero
        for i in range(m):
            c.A[i].pop(k, None)
        c.cost[k] = float(rng.randint(1, 3))  # min c x with c>0 over a free column: unbounded below along -x
    c.maximize = seed % 5 == 0
    if c.maximize:
        c.cost = [-x for x in c.cost]
    return Case("unbounded", seed, encode(c, rng, {}).render(), expect="unbounded")
