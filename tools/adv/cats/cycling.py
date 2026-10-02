"""cycling: Beale's and Kuhn's classic cycling examples (row/column permuted, power-of-two scaled) and homogeneous degenerate LPs that make textbook pivot rules cycle."""

from base import Case
from mpsio import Mps

NAME = "cycling"
COUNT = 100
SEED_BASE = 27000
MIP = False


def gen(rng, seed):
    v = seed % 3
    mp = Mps("CYC")
    mp.row("N", "OBJ")
    if v == 0:  # Beale (1955)
        spec = [("L", "R1", {"X4": 0.25, "X5": -8, "X6": -1, "X7": 9}),
                ("L", "R2", {"X4": 0.5, "X5": -12, "X6": -0.5, "X7": 3}),
                ("L", "R3", {"X6": 1})]
        cost = {"X4": -0.75, "X5": 20, "X6": -0.5, "X7": 6}
        rhs = {"R1": 0, "R2": 0, "R3": 1}
    elif v == 1:  # Kuhn's example with a box so the LP is bounded
        spec = [("L", "R1", {"X1": -2, "X2": -9, "X3": 1, "X4": 9}),
                ("L", "R2", {"X1": 1 / 3.0, "X2": 1, "X3": -1 / 3.0, "X4": -2})]
        cost = {"X1": -2, "X2": -3, "X3": 1, "X4": 12}
        rhs = {"R1": 0, "R2": 0}
    else:  # homogeneous rows (origin degenerate) with a normalising row
        nn = rng.randint(3, 8)
        mm = rng.randint(2, 6)
        spec = []
        for i in range(mm):
            spec.append(("L", "R%d" % i, {"X%d" % j: float(rng.randint(-4, 4)) for j in range(nn) if rng.random() < 0.7}))
        spec.append(("L", "RN", {"X%d" % j: 1.0 for j in range(nn)}))
        cost = {"X%d" % j: float(rng.randint(-5, 5)) for j in range(nn)}
        rhs = {r[1]: 0 for r in spec}
        rhs["RN"] = 1
    perm = list(range(len(spec)))
    rng.shuffle(perm)  # row order and column order matter to pivoting rules
    spec = [spec[p] for p in perm]
    scale = 2.0 ** rng.randint(-3, 3) if rng.random() < 0.5 else 1.0  # exact scaling
    for t, r, _ in spec:
        mp.row(t, r)
    cols = sorted(cost)
    rng.shuffle(cols)
    for c in cols:
        mp.add(c, "OBJ", cost[c] * scale)
        for t, r, d in spec:
            if c in d and d[c] != 0:
                mp.add(c, r, d[c])
    for r, b in rhs.items():
        mp.rhs.append((r, b))
    if v == 1:
        for c in cols:
            mp.bounds.append(("UP", c, 10))
    return Case("cycling", seed, mp.render())
