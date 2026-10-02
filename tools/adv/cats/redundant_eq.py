"""redundant_eq: redundant equalities (consistent, feasible) and inconsistent copies (infeasible)."""

from base import Case, dims
from core import rand_core, encode

NAME = "redundant_eq"
COUNT = 100
SEED_BASE = 34000
MIP = False


def gen(rng, seed):
    """Redundant equalities, sometimes with an inconsistent RHS (infeasible), sometimes duplicate with a free var."""
    m, n = dims(rng, 2, 6)
    c = rand_core(rng, m, n, rowkinds=("eq",), colkinds=("nonneg", "box", "free"))
    bad = seed % 3 == 0
    for _ in range(rng.randint(1, 3)):
        lam = [rng.randint(-2, 2) for _ in range(m)]
        if not any(lam):
            lam[0] = 1
        row = {}
        for l, i in zip(lam, range(m)):
            for j, v in c.A[i].items():
                row[j] = row.get(j, 0.0) + l * v
        row = {j: v for j, v in row.items() if v != 0}
        b = sum(l * c.row_lo[i] for l, i in zip(lam, range(m)))
        c.A.append(row)
        c.row_lo.append(b), c.row_up.append(b)
        c.m += 1
    expect = "optimal"
    if bad:
        c.row_lo[-1] += 1.0
        c.row_up[-1] += 1.0
        expect = "infeasible"
    return Case("redundant_eq", seed, encode(c, rng, {}).render(), expect=expect)
