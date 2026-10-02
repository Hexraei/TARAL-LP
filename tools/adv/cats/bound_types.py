"""bound_types: LO/UP/FX/FR/MI/PL combinations and orders on LP columns."""

from base import Case, dims
from core import rand_core, encode

NAME = "bound_types"
COUNT = 100
SEED_BASE = 42000
MIP = False


def gen(rng, seed):
    m, n = dims(rng, 3, 8)
    c = rand_core(rng, m, n, colkinds=("nonneg", "box", "free", "lowneg", "upneg", "fixed"))
    return Case("bound_types", seed, encode(c, rng, {}).render(), expect="optimal")
