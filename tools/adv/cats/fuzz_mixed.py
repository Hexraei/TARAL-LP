"""fuzz_mixed: random mixtures of the above (up to ~30x40), duplicates, zeros, ranges, constants, MAX."""

from base import Case, dims
from core import rand_core, encode

NAME = "fuzz_mixed"
COUNT = 250
SEED_BASE = 43000
MIP = False


def gen(rng, seed):
    m, n = dims(rng, 2, 12)
    c = rand_core(rng, m, n, dens=rng.choice([0.2, 0.5, 0.9]),
                  colkinds=("nonneg", "box", "free", "lowneg", "upneg", "fixed"),
                  rowkinds=("eq", "le", "ge", "rng"), slack0=rng.choice([0.2, 0.6, 1.0]))
    c.const = float(rng.choice([0, 0, 5, -2.5]))
    c.maximize = rng.random() < 0.3
    if c.maximize:
        c.cost = [-x for x in c.cost]
    st = {"split": rng.choice([0, 0, 0.3]), "zeros": rng.choice([0, 0, 2]),
          "sense": rng.choice(["section", "inline"])}
    if st["split"]:
        st["coef"] = "int"
    return Case("fuzz_mixed", seed, encode(c, rng, st).render(), expect="optimal")
