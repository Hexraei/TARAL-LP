"""Shared case container and size helper for the category generators."""


class Case:
    """One generated MPS case. `expect` is set only where the construction proves the status; `exact_obj` only
    where the optimum is known exactly (used to arbitrate when the HiGHS reference is itself inexact)."""

    def __init__(self, cat, seed, text, mip=False, expect=None, exact_obj=None, note=""):
        self.cat, self.seed, self.text, self.mip = cat, seed, text, mip
        self.expect, self.exact_obj, self.note = expect, exact_obj, note
        self.id = "%s-%d" % (cat, seed)


def dims(rng, lo=2, hi=9):
    """Mostly small models; one in four is 'large' (up to ~30 x 40) so ties/pivot sequences get long."""
    if rng.random() < 0.25:
        return rng.randint(10, 30), rng.randint(15, 40)
    m = rng.randint(lo, hi)
    n = rng.randint(lo, hi + 3)
    return m, n
