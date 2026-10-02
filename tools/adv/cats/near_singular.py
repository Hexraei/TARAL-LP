"""near_singular: square equality systems with condition number 1e8-1e12 (prescribed singular values, Hilbert, nearly dependent rows), plain and objective-stable costs; an exact rational optimum arbitrates."""
import random
from fractions import Fraction

import numpy as np

from base import Case
from mpsio import Mps

NAME = "near_singular"
COUNT = 150
SEED_BASE = 37000
MIP = False


def _exact_solve(B, b):
    m = len(B)
    M = [[Fraction(x) for x in row] + [Fraction(bb)] for row, bb in zip(B, b)]
    for c in range(m):
        p = next((r for r in range(c, m) if M[r][c] != 0), None)
        if p is None:
            return None
        M[c], M[p] = M[p], M[c]
        for r in range(m):
            if r != c and M[r][c] != 0:
                f = M[r][c] / M[c][c]
                M[r] = [a - f * bb for a, bb in zip(M[r], M[c])]
    return [M[i][m] / M[i][i] for i in range(m)]


def gen(rng, seed):
    v = seed % 3
    rs = np.random.RandomState(seed)
    m = rng.randint(3, 8) if v != 1 else rng.randint(4, 8)
    if v == 0:  # prescribed singular values: cond in [1e8, 1e12]
        kappa = 10.0 ** rng.uniform(8, 12)
        U, _ = np.linalg.qr(rs.randn(m, m))
        V, _ = np.linalg.qr(rs.randn(m, m))
        s = np.logspace(0, -np.log10(kappa), m)
        B = (U * s) @ V.T
        B = [[float(x) for x in row] for row in B]
    elif v == 1:  # Hilbert matrix (cond 1e7..1e11 for m = 6..9)
        m = rng.randint(6, 9)
        B = [[1.0 / (i + j + 1) for j in range(m)] for i in range(m)]
    else:  # nearly dependent rows: row_k = row_0 + eps * noise
        eps = 10.0 ** -rng.randint(8, 12)
        B = [[float(rng.randint(-5, 5)) or 1.0 for _ in range(m)] for _ in range(m)]
        B[-1] = [B[0][j] + eps * rng.uniform(-1, 1) for j in range(m)]
    x0 = [float(rng.randint(1, 4)) for _ in range(m)]
    b = [float(sum(Fraction(B[i][j]) * Fraction(x0[j]) for j in range(m))) for i in range(m)]
    cost = [float(rng.randint(1, 5)) for _ in range(m)]
    if (seed // 3) % 2:  # objective-stable variant: c = B'y, so c'x = y'b whatever error x carries along the near-null space
        y = [float(rng.randint(1, 3)) for _ in range(m)]
        cost = [float(sum(Fraction(B[i][j]) * Fraction(y[i]) for i in range(m))) for j in range(m)]
    mp = Mps("NSING")
    mp.row("N", "OBJ")
    for i in range(m):
        mp.row("E", "R%d" % i)
    for j in range(m):
        mp.add("X%d" % j, "OBJ", cost[j])
        for i in range(m):
            mp.add("X%d" % j, "R%d" % i, B[i][j])
    for i in range(m):
        mp.rhs.append(("R%d" % i, b[i]))
    ex = _exact_solve(B, b)
    exact = None
    note = ""
    if ex is not None and all(t >= 0 for t in ex):
        exact = float(sum(Fraction(cost[j]) * ex[j] for j in range(m)))
        note = "square system: unique feasible point, exact objective known"
    return Case("near_singular", seed, mp.render(), exact_obj=exact, expect="optimal" if exact is not None else None,
                note=note)
