"""A feasible-by-construction random LP/MILP 'core' and its MPS encodings, shared by all generators.

The core stores the *intended* model (row ranges, column bounds, costs, integrality). encode() writes
it as MPS choosing among the equivalent spellings the format allows (row type + RHS + RANGES, bound
types and order, duplicate coefficient lines, infinity spellings, ...). The oracle parser re-derives the
model from the text; `selftest` in gen.py checks that both agree, which validates writer and oracle.
"""
from mpsio import INF, Mps, Model


class Core:
    def __init__(self, m, n):
        self.m, self.n = m, n
        self.A = [dict() for _ in range(m)]  # row -> {col: value}
        self.row_lo = [-INF] * m
        self.row_up = [INF] * m
        self.col_lo = [0.0] * n
        self.col_up = [INF] * n
        self.cost = [0.0] * n
        self.const = 0.0
        self.maximize = False
        self.ints = set()
        self.x0 = [0.0] * n  # a feasible point (when the core is feasible)

    def model(self):
        """Intended model, in the oracle's own data structure (for the writer/oracle self-test)."""
        md = Model()
        md.maximize = self.maximize
        md.obj_const = self.const
        md.row_names = ["R%d" % i for i in range(self.m)]
        md.row_lo, md.row_up = list(self.row_lo), list(self.row_up)
        md.col_names = ["C%d" % j for j in range(self.n)]
        md.cost = list(self.cost)
        md.col_lo, md.col_up = list(self.col_lo), list(self.col_up)
        md.is_int = [j in self.ints for j in range(self.n)]
        md.entries = [dict() for _ in range(self.n)]
        for i in range(self.m):
            for j, v in self.A[i].items():
                md.entries[j][i] = v
        return md


def rand_core(rng, m, n, dens=0.4, coef="int", colkinds=("nonneg", "box"), rowkinds=("eq", "le", "ge", "rng"),
              slack0=0.4, cost="dual", integer=None, xrange=5):
    """Random core with a known feasible point x0. cost='dual' builds c = d + A'y so that the LP is bounded."""
    c = Core(m, n)
    for j in range(n):
        k = rng.choice(colkinds)
        if k == "nonneg":
            lo, up = 0.0, INF
        elif k == "box":
            lo = float(rng.choice([0, 0, rng.randint(-3, 3)]))
            up = lo + rng.randint(1, 10)
        elif k == "free":
            lo, up = -INF, INF
        elif k == "lowneg":
            lo, up = float(rng.randint(-5, 0)), INF
        elif k == "upneg":
            lo, up = -INF, float(rng.randint(-3, 5))
        elif k == "fixed":
            lo = up = float(rng.randint(-3, 5))
        else:
            raise ValueError(k)
        c.col_lo[j], c.col_up[j] = lo, up
        a = lo if lo > -INF else (up - rng.randint(0, xrange) if up < INF else -rng.randint(0, xrange))
        b = up if up < INF else a + rng.randint(0, xrange)
        if rng.random() < 0.35:
            x = a if rng.random() < 0.5 else b  # on a bound: primal degeneracy
        else:
            x = float(rng.randint(int(a), int(b)))
        c.x0[j] = float(x)
    if integer:
        c.ints = {j for j in range(n) if rng.random() < integer}
        for j in c.ints:  # integer columns keep integer bounds (all generated bounds already are)
            pass

    def rcoef():
        if coef == "dyadic":
            v = rng.randint(-24, 24) / 8.0
            return v if v else 0.125
        v = rng.randint(-5, 5)
        return float(v if v else rng.choice([-1, 1]))

    for i in range(m):
        cols = [j for j in range(n) if rng.random() < dens]
        if not cols:
            cols = [rng.randrange(n)]
        for j in cols:
            c.A[i][j] = rcoef()
    for i in range(m):
        act = sum(v * c.x0[j] for j, v in c.A[i].items())
        k = rng.choice(rowkinds)
        s1 = 0.0 if rng.random() < slack0 else float(rng.randint(1, 5))
        s2 = 0.0 if rng.random() < slack0 else float(rng.randint(1, 5))
        if k == "eq":
            c.row_lo[i] = c.row_up[i] = act
        elif k == "le":
            c.row_up[i] = act + s1
        elif k == "ge":
            c.row_lo[i] = act - s1
        else:
            c.row_lo[i], c.row_up[i] = act - s1, act + s2
    if cost == "dual":
        y = [0.0] * m
        for i in range(m):
            lo, up = c.row_lo[i], c.row_up[i]
            if rng.random() < 0.6:
                mag = float(rng.randint(1, 3))
                if lo > -INF and up < INF:
                    y[i] = mag * rng.choice([-1, 1])
                elif up < INF:
                    y[i] = -mag  # L row: y <= 0 for min
                else:
                    y[i] = mag
        cc = [0.0] * n
        for i in range(m):
            for j, v in c.A[i].items():
                cc[j] += v * y[i]
        for j in range(n):
            lo, up = c.col_lo[j], c.col_up[j]
            if lo == up:
                d = float(rng.randint(-3, 3))
            elif lo > -INF and up < INF:
                d = float(rng.randint(-3, 3))
            elif lo > -INF:
                d = float(rng.choice([0, 0, 1, 2, 3]))
            elif up < INF:
                d = -float(rng.choice([0, 0, 1, 2, 3]))
            else:
                d = 0.0
            c.cost[j] = cc[j] + d
    else:
        c.cost = [float(rng.randint(-5, 5)) for _ in range(n)]
    return c


# --------------------------------------------------------------------------------------------- encoding
def _row_encodings(lo, up, rng):
    """All (type, rhs, range-or-None) spellings of [lo, up]; returns one chosen at random."""
    opts = []
    if lo == up:
        opts.append(("E", lo, None))
        opts.append(("L", up, 0.0))
        opts.append(("G", lo, 0.0))
        opts.append(("E", lo, 0.0))
        opts.append(("E", lo, "-0"))
    elif lo == -INF and up == INF:
        return ("N", 0.0, None)
    elif lo == -INF:
        opts.append(("L", up, None))
    elif up == INF:
        opts.append(("G", lo, None))
    else:
        w = up - lo
        opts += [("L", up, w), ("G", lo, w), ("E", lo, w), ("E", up, -w)]
        opts += [("L", up, -w), ("G", lo, -w)]  # sign of RANGES on L/G rows is ignored (|R|)
    return rng.choice(opts)


def _special(core_lo, core_up, style):
    """True when an integer column is marked integer by its BV/LI/UI bound lines instead of a MARKER block."""
    return style.get("int_bounds") == "special" and not (core_lo == -INF and core_up == INF)


def _emit_bounds(mp, col, lo, up, is_int, rng, style):
    B = mp.bounds
    huge = style.get("huge")  # infinity spelling for bounds
    posinf = {"1e20": "1e20", "1e30": "1e30", "inf": "1e+30"}[huge] if huge else None
    if is_int and _special(lo, up, style):
        if lo == 0 and up == 1:
            B.append(("BV", col, None))
            return
        B.append(("LI", col, lo) if lo > -INF else ("MI", col, None))
        if up < INF:
            B.append(("UI", col, up))
        return
    if lo == -INF and up == INF:
        r = rng.random()
        if huge and r < 0.5:
            B.append(("LO", col, "-" + posinf))
            B.append(("UP", col, posinf))
        elif r < 0.8:
            B.append(("FR", col, None))
        else:
            B.append(("MI", col, None))
            B.append(("PL", col, None))
        return
    if lo == up:
        if rng.random() < 0.7 or is_int:
            B.append(("FX", col, lo))
        else:
            B.append(("LO", col, lo))
            B.append(("UP", col, up))
        return
    if lo == -INF:
        B.append(("MI", col, None) if not huge or rng.random() < 0.5 else ("LO", col, "-" + posinf))
        B.append(("UP", col, up))
        return
    if lo != 0 or rng.random() < 0.1:
        B.append(("LO", col, lo))
    if up != INF:
        B.append(("UP", col, up))
    elif huge and rng.random() < 0.5:
        B.append(("UP", col, posinf))
    elif is_int or rng.random() < 0.05:
        B.append(("PL", col, None))


def encode(core, rng, style=None):
    """Return an Mps for the core. style keys: split (prob. of splitting an entry in two lines), zeros (extra explicit
    0 entries), negzero, huge (None|'1e20'|'1e30'|'inf'), sense ('section'|'inline'), int_bounds ('marker'|'special'),
    rhs_const_row (bool), extra_free_rows."""
    st = style or {}
    mp = Mps("ADV")
    mp.row("N", "OBJ")
    names = []
    for i in range(core.m):
        t, rhs, R = _row_encodings(core.row_lo[i], core.row_up[i], rng)
        nm = "R%d" % i
        mp.row(t, nm)
        names.append(nm)
        if t != "N":
            if rhs != 0 or st.get("zero_rhs_lines") or rng.random() < 0.2:
                mp.rhs.append((nm, rhs if not (st.get("negzero") and rhs == 0) else "-0"))
        if R is not None and t != "N":
            mp.ranges.append((nm, R))
    huge = st.get("huge")
    if huge and rng.random() < 1.0:  # extra rows whose RHS is the infinite sentinel
        pass
    for j in range(core.n):
        cn = "C%d" % j
        mp.col(cn, integer=(j in core.ints and not _special(core.col_lo[j], core.col_up[j], st)))
    ent = {j: [] for j in range(core.n)}
    for j in range(core.n):
        if core.cost[j] != 0 or st.get("zero_cost_lines") or rng.random() < 0.1:
            ent[j].append(("OBJ", core.cost[j] if not (st.get("negzero") and core.cost[j] == 0) else "-0"))
    for i in range(core.m):
        for j, v in core.A[i].items():
            sp = st.get("split", 0.0)
            if rng.random() < sp:
                if st.get("coef") == "dyadic" or float(v).is_integer():
                    a = float(rng.randint(-8, 8))
                    if st.get("coef") == "dyadic":
                        a = rng.randint(-16, 16) / 8.0
                    ent[j].append((names[i], a))
                    ent[j].append((names[i], v - a))
                else:
                    ent[j].append((names[i], v))
            else:
                ent[j].append((names[i], v))
    nz = st.get("zeros", 0)
    for _ in range(nz):  # explicit zeros on rows (and the objective)
        j = rng.randrange(core.n)
        r = rng.choice(names + ["OBJ"]) if names else "OBJ"
        ent[j].append((r, "-0" if st.get("negzero") else 0))
    for j in range(core.n):
        e = ent[j]
        if st.get("shuffle_entries", True):
            rng.shuffle(e)
        if not e:
            e = [("OBJ", 0)]
        for r, v in e:
            mp.add("C%d" % j, r, v)
    if core.const != 0 or st.get("zero_const_line"):
        mp.rhs.insert(0, ("OBJ", -core.const if core.const else ("-0" if st.get("negzero") else 0)))
    if core.maximize:
        mp.objsense = "MAX"
        mp.objsense_style = st.get("sense", "section")
    elif st.get("explicit_min"):
        mp.objsense = "MIN"
        mp.objsense_style = st.get("sense", "section")
    for j in range(core.n):
        _emit_bounds(mp, "C%d" % j, core.col_lo[j], core.col_up[j], j in core.ints, rng, st)
    return mp
