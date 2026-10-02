"""neg_lower: negative lower bounds, negative upper bounds (MI+UP), boxes straddling zero, free and fixed columns."""

from base import Case, dims
from core import rand_core, encode

NAME = "neg_lower"
COUNT = 100
SEED_BASE = 41000
MIP = False


def gen(rng, seed):
    m, n = dims(rng)
    c = rand_core(rng, m, n, colkinds=("lowneg", "upneg", "box", "free", "fixed"))
    return Case("neg_lower", seed, encode(c, rng, {}).render(), expect="optimal")
