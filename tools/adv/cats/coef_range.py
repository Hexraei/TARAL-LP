"""coef_range: exactly rescaled models (row scales 10^-4..10^4, column scales 10^-5..10^5) so entries span 1e-9..1e9 while the optimum is unchanged."""

from base import Case, dims
from core import Core, rand_core, encode
from mpsio import INF

NAME = "coef_range"
COUNT = 200
SEED_BASE = 38000
MIP = False


def gen(rng, seed):
    """Row scales 10^r_i and column scales 10^s_j make the model exactly equivalent (x'_j = x_j / 10^s_j);
    entries span 1e-9..1e9. Half the cases sit at the extremes of the exponent range."""
    m, n = dims(rng)
    c = rand_core(rng, m, n, colkinds=("nonneg", "box", "lowneg"))
    if seed % 2:
        r = [rng.choice([-4, 4, 0]) for _ in range(m)]
        s = [rng.choice([-5, 5, 0]) for _ in range(n)]
    else:
        r = [rng.randint(-4, 4) for _ in range(m)]
        s = [rng.randint(-5, 5) for _ in range(n)]
    sc = Core(m, n)
    for i in range(m):
        for j, v in c.A[i].items():
            sc.A[i][j] = v * 10.0 ** (r[i] + s[j])
        f = 10.0 ** r[i]
        sc.row_lo[i] = c.row_lo[i] * f if c.row_lo[i] > -INF else -INF
        sc.row_up[i] = c.row_up[i] * f if c.row_up[i] < INF else INF
    for j in range(n):
        g = 10.0 ** s[j]
        sc.col_lo[j] = c.col_lo[j] / g if c.col_lo[j] > -INF else -INF
        sc.col_up[j] = c.col_up[j] / g if c.col_up[j] < INF else INF
        sc.cost[j] = c.cost[j] * g
        sc.x0[j] = c.x0[j] / g
    return Case("coef_range", seed, encode(sc, rng, {}).render(), expect="optimal")
