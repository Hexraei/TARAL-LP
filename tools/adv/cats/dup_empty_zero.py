"""dup_empty_zero: duplicate / scaled / sign-flipped rows, empty rows (satisfied and violated), explicit zero entries, zero rows, empty and zero columns, duplicate columns."""

from base import Case, dims
from core import rand_core, encode
from mpsio import INF

NAME = "dup_empty_zero"
COUNT = 200
SEED_BASE = 29000
MIP = False


def gen(rng, seed):
    m, n = dims(rng, 2, 7)
    c = rand_core(rng, m, n)
    v = seed % 6
    expect = "optimal"
    st = {}
    if v == 0:  # exact duplicate rows (and scaled / sign-flipped copies)
        for _ in range(rng.randint(1, 3)):
            i = rng.randrange(m)
            k = rng.choice([1.0, 2.0, -1.0, 0.5, -4.0])
            c.A.append({j: a * k for j, a in c.A[i].items()})
            lo, up = c.row_lo[i], c.row_up[i]
            nl, nu = (lo * k, up * k) if k > 0 else (up * k, lo * k)
            c.row_lo.append(nl)
            c.row_up.append(nu)
            c.m += 1
    elif v == 1:  # empty rows, satisfied
        for _ in range(rng.randint(1, 3)):
            c.A.append({})
            kind = rng.choice(["E", "L", "G"])
            if kind == "E":
                c.row_lo.append(0.0), c.row_up.append(0.0)
            elif kind == "L":
                c.row_lo.append(-INF), c.row_up.append(float(rng.randint(0, 3)))
            else:
                c.row_lo.append(-float(rng.randint(0, 3))), c.row_up.append(INF)
            c.m += 1
    elif v == 2:  # empty row, violated -> infeasible
        c.A.append({})
        kind = rng.choice(["E", "L", "G"])
        if kind == "E":
            b = float(rng.choice([1, -1, 3]))
            c.row_lo.append(b), c.row_up.append(b)
        elif kind == "L":
            c.row_lo.append(-INF), c.row_up.append(-float(rng.randint(1, 3)))
        else:
            c.row_lo.append(float(rng.randint(1, 3))), c.row_up.append(INF)
        c.m += 1
        expect = "infeasible"
    elif v == 3:  # explicit zero entries, zero rows, zero columns
        st["zeros"] = rng.randint(2, 8)
        c.A.append({})
        c.row_lo.append(-INF), c.row_up.append(float(rng.randint(0, 2)))
        c.m += 1
    elif v == 4:  # columns without entries: cost >= 0 harmless, cost < 0 with no upper bound unbounded
        k = rng.randint(1, 3)
        neg = rng.random() < 0.5
        for _ in range(k):
            c.n += 1
            c.col_lo.append(0.0), c.col_up.append(INF), c.x0.append(0.0)
            c.cost.append(-float(rng.randint(1, 3)) if neg else float(rng.randint(0, 3)))
        expect = "unbounded" if neg else "optimal"
    else:  # duplicate identical columns (same entries) and a column that is a multiple of another
        j = rng.randrange(n)
        for k in (1.0, 1.0, 2.0):
            c.n += 1
            c.col_lo.append(min(c.col_lo[j], 0.0)), c.col_up.append(max(c.col_up[j], 0.0))
            c.x0.append(0.0)
            c.cost.append(c.cost[j] * k)
            for i in range(c.m):
                if j in c.A[i]:
                    c.A[i][c.n - 1] = c.A[i][j] * k
    txt = encode(c, rng, st).render()
    return Case("dup_empty_zero", seed, txt, expect=expect)
