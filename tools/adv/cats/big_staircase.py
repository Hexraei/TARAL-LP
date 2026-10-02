"""big_staircase: multi-period production planning (staircase structure): T periods x K products, inventory balance
equalities linking consecutive periods, shared capacity rows, 500-3000 rows. Capacity covers demand by construction
(optimal); a few instances have overtime columns with zero cost (dual degeneracy) and ranges on capacity rows."""

from base import Case
from mpsio import Mps

NAME = "big_staircase"
COUNT = 24
SEED_BASE = 61000
MIP = False
BIG = True


def gen(rng, seed):
    K = rng.randint(3, 12)
    T = max(5, int(round(500 * (6.0 ** rng.random()))) // (K + 1))
    mp = Mps("STAIR")
    mp.row("N", "OBJ")
    dem = [[rng.randint(0, 9) for _ in range(K)] for _ in range(T)]
    for t in range(T):
        for k in range(K):
            mp.row("E", "B%d_%d" % (t, k))  # I[t-1] + P[t] - I[t] = demand
        mp.row("L", "C%d" % t)  # sum of production hours <= capacity (+ overtime)
    use = [rng.randint(1, 3) for _ in range(K)]
    for t in range(T):
        cap = sum(dem[t][k] * use[k] for k in range(K)) + rng.randint(0, 6) * rng.choice([0, 1, 1, 2])
        mp.rhs.append(("C%d" % t, cap + sum(use) * rng.choice([0, 1, 2])))  # slack for catching up through stock
        if rng.random() < 0.3:
            mp.ranges.append(("C%d" % t, rng.randint(1, 20)))
        for k in range(K):
            p, i = "P%d_%d" % (t, k), "I%d_%d" % (t, k)
            mp.add(p, "OBJ", rng.choice([1, 2, 2, 3]))
            mp.add(p, "B%d_%d" % (t, k), 1)
            mp.add(p, "C%d" % t, use[k])
            mp.add(i, "OBJ", rng.choice([0, 1, 1, 2]))
            mp.add(i, "B%d_%d" % (t, k), -1)
            if t + 1 < T:
                mp.add(i, "B%d_%d" % (t + 1, k), 1)
            mp.rhs.append(("B%d_%d" % (t, k), dem[t][k]))
        if rng.random() < 0.25:  # free overtime column (zero cost), bounded by a box
            o = "O%d" % t
            mp.add(o, "C%d" % t, -1)
            mp.bounds.append(("UP", o, rng.randint(1, 6)))
    return Case("big_staircase", seed, mp.render(), expect="optimal", note="T=%d K=%d" % (T, K))
