"""rank_deficient: equality systems whose rank is below the row count: consistent dependent rows (feasible)."""

from base import Case
from core import Core, encode

NAME = "rank_deficient"
COUNT = 150
SEED_BASE = 33000
MIP = False


def gen(rng, seed):
    """Equalities with rank below the row count: dependent rows with consistent RHS (feasible)."""
    k = rng.randint(2, 5)
    n = rng.randint(k + 1, 9)
    x0 = [float(rng.randint(0, 4)) for _ in range(n)]
    base = [{j: float(rng.randint(-3, 3)) for j in range(n)} for _ in range(k)]
    rows, rhs = [], []
    for b in base:
        rows.append(b)
        rhs.append(sum(v * x0[j] for j, v in b.items()))
    for _ in range(rng.randint(1, 4)):
        lam = [rng.randint(-2, 3) for _ in base]
        rows.append({j: sum(l * b[j] for l, b in zip(lam, base)) for j in range(n)})
        rhs.append(sum(l * r for l, r in zip(lam, rhs[:k])))
    order = list(range(len(rows)))
    rng.shuffle(order)
    c = Core(len(rows), n)
    c.x0 = x0
    for t, p in enumerate(order):
        c.A[t] = {j: v for j, v in rows[p].items() if v != 0}
        c.row_lo[t] = c.row_up[t] = rhs[p]
    c.col_up = [x0[j] + rng.randint(1, 5) for j in range(n)]  # boxed -> bounded
    c.cost = [float(rng.randint(-4, 4)) for _ in range(n)]
    return Case("rank_deficient", seed, encode(c, rng, {}).render(), expect="optimal")
