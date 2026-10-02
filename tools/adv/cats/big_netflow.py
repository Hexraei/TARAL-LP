"""big_netflow: min-cost network flow on random sparse digraphs with 500-3000 nodes (rows) and about 4 arcs per node.
Node balances are equalities whose rows sum to zero, so the system has rank (nodes - 1): one redundant equality, plus
heavy cost ties (integer costs 1..5) and capacities that make many arcs tight (primal and dual degeneracy). A
high-capacity bidirectional backbone ring keeps it feasible: expected status is optimal."""

from base import Case
from mpsio import Mps

NAME = "big_netflow"
COUNT = 24
SEED_BASE = 60000
MIP = False
BIG = True


def gen(rng, seed):
    nv = int(round(500 * (6.0 ** rng.random())))
    arcs = []
    perm = list(range(nv))
    rng.shuffle(perm)
    tot = 0
    sup = [0] * nv
    for _ in range(nv // 4):  # balanced supplies: sources and sinks of 1..9 units
        a, b = rng.randrange(nv), rng.randrange(nv)
        q = rng.randint(1, 9)
        sup[a] += q
        sup[b] -= q
        tot += q
    for k in range(nv):  # backbone ring, both directions, capacity above the total supply
        a, b = perm[k], perm[(k + 1) % nv]
        arcs.append((a, b, rng.randint(1, 5), tot + 1))
        arcs.append((b, a, rng.randint(1, 5), tot + 1))
    for _ in range(2 * nv):
        a, b = rng.randrange(nv), rng.randrange(nv)
        if a != b:
            arcs.append((a, b, rng.choice([1, 1, 2, 2, 3, 4, 5]), rng.choice([1, 2, 3, 5, 9, 20])))
    mp = Mps("NETFLOW")
    mp.row("N", "OBJ")
    for v in range(nv):
        mp.row("E", "N%d" % v)
    for k, (a, b, cost, cap) in enumerate(arcs):
        cn = "F%d" % k
        mp.add(cn, "OBJ", cost)
        mp.add(cn, "N%d" % a, 1)  # leaves a
        mp.add(cn, "N%d" % b, -1)  # enters b
        mp.bounds.append(("UP", cn, cap))
    for v in range(nv):
        if sup[v]:
            mp.rhs.append(("N%d" % v, sup[v]))
    return Case("big_netflow", seed, mp.render(), expect="optimal", note="%d nodes, %d arcs" % (nv, len(arcs)))
