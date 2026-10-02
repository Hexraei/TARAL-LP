"""Seeded adversarial case registry. gen(category, i) -> Case, deterministic.

Each module in tools/adv/cats/ defines NAME, COUNT, SEED_BASE, MIP and gen(rng, seed). Case i of a category uses
seed SEED_BASE + i (SEED_BASE = 26000 + 1000*k, fixed per category so adding a category never moves existing seeds).
Only the stdlib random.Random(seed) (and numpy's legacy RandomState for dense linear algebra) is used.
`python3 tools/adv/gen.py --seeds` prints the seed list committed as tests/seeds.txt.
"""
import importlib
import os
import pkgutil
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cats

_MODS = sorted((importlib.import_module("cats." + n) for _, n, _ in pkgutil.iter_modules(cats.__path__)),
               key=lambda m: m.SEED_BASE)
# Modules with BIG = True (models of hundreds to a few thousand rows) are a separate group: they are not part of the
# base suite (ALL / LP_CATS / MILP_CATS, tests/seeds.txt) and are selected with harness.py --kind big.
_BASE = [m for m in _MODS if not getattr(m, "BIG", False)]
_BIG = [m for m in _MODS if getattr(m, "BIG", False)]
LP_CATS = [(m.NAME, m.gen, m.COUNT) for m in _BASE if not m.MIP]
MILP_CATS = [(m.NAME, m.gen, m.COUNT) for m in _BASE if m.MIP]
ALL = [(m.NAME, m.gen, m.COUNT) for m in _BASE]
BIG_CATS = [(m.NAME, m.gen, m.COUNT) for m in _BIG]
BASE = {m.NAME: m.SEED_BASE for m in _MODS}
GEN = {m.NAME: m.gen for m in _MODS}
COUNT = {m.NAME: m.COUNT for m in _MODS}
IS_MIP = {m.NAME: True for m in _MODS if m.MIP}


def gen(cat, i):
    seed = BASE[cat] + i
    case = GEN[cat](random.Random(seed), seed)
    case.mip = case.mip or IS_MIP.get(cat, False)
    return case


if __name__ == "__main__":
    if "--seeds" in sys.argv or "--seeds-big" in sys.argv:
        print("# category  kind  first_seed  last_seed  count   (case i has seed first_seed + i; fixed, deterministic)")
        for m in (_BIG if "--seeds-big" in sys.argv else _BASE):
            print("%-16s %-5s %d %d %d" % (m.NAME, "MILP" if m.MIP else "LP", m.SEED_BASE, m.SEED_BASE + m.COUNT - 1, m.COUNT))
        grp = _BIG if "--seeds-big" in sys.argv else _BASE
        lp = sum(m.COUNT for m in grp if not m.MIP)
        mi = sum(m.COUNT for m in grp if m.MIP)
        print("# total LP %d  MILP %d" % (lp, mi))
