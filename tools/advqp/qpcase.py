"""Convex / nonconvex QP adversarial cases: data model, MPS writer and generators.

Model of a case (the generator's own numbers are the oracle; no engine or HiGHS reader is involved in defining it):
    sense  c'x + 0.5 x'Qx + const     (sense = min, or max when maximize)
    rlo <= A x <= rup,  lo <= x <= up       (|bound| >= 1e20 is infinite)
Q is the full symmetric matrix. A max problem with concave Q is convex.

Construction knowledge a generator may attach (never taken from a solver):
    expect_status   'optimal' | 'infeasible' | 'unbounded' | None (None = ask the reference)
    x0, expect_obj  a point proved optimal by an explicit KKT construction (see kkt_cost)
    ray             recession direction for 'unbounded' cases (integer data, exactly verifiable)
Seeds: category k uses base SEED_BASE + 1000*k and case i uses base + i (numpy RandomState only).
"""
import math

import numpy as np

INF = math.inf
SEED_BASE = 70000


# ------------------------------------------------------------------------------------------------ data model
class QCase:
    def __init__(self, cat, seed, Q, c, A, rlo, rup, lo, up, maximize=False, const=0.0, qform="lower", rowstyle=0):
        self.cat, self.seed = cat, seed
        self.id = "%s-%d" % (cat, seed)
        self.Q = np.asarray(Q, float)
        self.c = np.asarray(c, float)
        n = len(self.c)
        self.A = np.asarray(A, float).reshape(-1, n)
        self.rlo, self.rup = np.asarray(rlo, float), np.asarray(rup, float)
        self.lo, self.up = np.asarray(lo, float), np.asarray(up, float)
        self.maximize, self.const, self.qform, self.rowstyle = maximize, float(const), qform, rowstyle
        self.expect_status = None
        self.x0 = None
        self.expect_obj = None
        self.ray = None
        self.note = ""
        self.qdups = False   # split every Q entry into two summing lines
        self.qzeros = False  # insert explicit zero entries / cancelling pairs
        self.qempty = False  # a QUADOBJ header with no lines when Q is zero
        for a in ("rlo", "rup", "lo", "up"):  # 1e20 and beyond is infinity, in the data model too
            v = getattr(self, a)
            v = np.where(v >= 1e20, INF, np.where(v <= -1e20, -INF, v))
            setattr(self, a, v)

    @property
    def n(self):
        return len(self.c)

    @property
    def m(self):
        return self.A.shape[0]

    def obj(self, x):
        x = np.asarray(x, float)
        return float(self.c @ x + 0.5 * x @ self.Q @ x + self.const)

    def qmin(self):
        return -self.Q if self.maximize else self.Q

    def eig(self):
        """(lam_min, lam_max_abs) of the Hessian in the minimisation sense."""
        if self.n == 0:
            return 0.0, 0.0
        w = np.linalg.eigvalsh(self.qmin())
        return float(w[0]), float(np.max(np.abs(w)))


# ------------------------------------------------------------------------------------------------ MPS writer
def num(v):
    v = float(v)
    return str(int(v)) if v == int(v) and abs(v) < 1e15 else repr(v)


def write_mps(case, rs=None):
    """MPS text for the case. Style variants (row type spelling, bound spelling, Q listing) come from case fields."""
    n, m = case.n, case.m
    rs = rs or np.random.RandomState(case.seed)
    if np.any(case.rlo > case.rup):
        raise ValueError("empty row range is not representable in MPS (RANGES uses |R| for L/G rows)")
    out = ["NAME %s" % case.id.replace("-", "_")]
    if case.maximize:
        out += ["OBJSENSE", "    MAX"]
    out.append("ROWS")
    out.append(" N OBJ")
    kinds, rhs, rng = [], [], {}
    for i in range(m):
        lo, up = case.rlo[i], case.rup[i]
        if lo == up:
            k, b = "E", lo
        elif lo == -INF and up == INF:
            k, b = "N", 0.0  # free row: written as an extra N row (dropped by readers)
        elif lo == -INF:
            k, b = "L", up
        elif up == INF:
            k, b = "G", lo
        else:
            s = (case.rowstyle + i) % 3
            if s == 0:
                k, b, rng[i] = "L", up, up - lo
            elif s == 1:
                k, b, rng[i] = "G", lo, up - lo
            else:
                k, b, rng[i] = "E", lo, up - lo  # E with R>0 is [rhs, rhs+R]
        kinds.append(k)
        rhs.append(b)
        out.append(" %s R%d" % (k, i))
    out.append("COLUMNS")
    for j in range(n):
        out.append("    C%d OBJ %s" % (j, num(case.c[j])))
        for i in range(m):
            if case.A[i, j] != 0:
                out.append("    C%d R%d %s" % (j, i, num(case.A[i, j])))
    out.append("RHS")
    if case.const != 0:
        out.append("    RHS OBJ %s" % num(-case.const))  # objective constant = - RHS of the objective row
    for i in range(m):
        if kinds[i] != "N" and rhs[i] != 0:
            out.append("    RHS R%d %s" % (i, num(rhs[i])))
    if rng:
        out.append("RANGES")
        for i, r in rng.items():
            out.append("    RNG R%d %s" % (i, num(r)))
    out.append("BOUNDS")
    for j in range(n):
        lo, up = case.lo[j], case.up[j]
        if lo == up:
            out.append(" FX BND C%d %s" % (j, num(lo)))
        elif lo == -INF and up == INF:
            out.append(" FR BND C%d" % j)
        else:
            if lo == -INF:
                out.append(" MI BND C%d" % j)
            if up != INF:  # UP before LO: UP with a negative value and a default lower bound sets lo=-inf first
                out.append(" UP BND C%d %s" % (j, num(up)))
            if lo != -INF and (lo != 0 or up < 0):
                out.append(" LO BND C%d %s" % (j, num(lo)))
    Q = case.Q
    ent = [(i, j, Q[i, j]) for j in range(n) for i in range(j, n) if Q[i, j] != 0]
    if case.qform == "full":
        out.append("QMATRIX")
        for i, j, v in ent:
            out.append("    C%d C%d %s" % (j, i, num(v)))
            if i != j:
                out.append("    C%d C%d %s" % (i, j, num(v)))
    elif case.qform == "qsection":
        out.append("QSECTION OBJ")
        for i, j, v in ent:
            out.append("    C%d C%d %s" % (j, i, num(v)))
            if i != j:
                out.append("    C%d C%d %s" % (i, j, num(v)))
    else:
        if ent or case.qempty or case.qzeros:
            out.append("QUADOBJ")
        lines = []
        for i, j, v in ent:
            a, b = (j, i) if case.qform == "lower" else (i, j)  # token order only; one triangle per entry
            if case.qdups:
                p = float(np.round(rs.uniform(-2, 2), 3))
                lines.append((a, b, p))
                lines.append((a, b, v - p))
            else:
                lines.append((a, b, v))
        if case.qzeros:
            for _ in range(max(1, n // 2)):
                a, b = int(rs.randint(n)), int(rs.randint(n))
                lines.append((a, b, 0.0))
            if ent:
                i, j, v = ent[int(rs.randint(len(ent)))]
                t = float(rs.randint(1, 9))
                lines.append((j, i, t))  # cancelling pair on one entry
                lines.append((j, i, -t))
        if case.qdups or case.qzeros:
            order = rs.permutation(len(lines))
            lines = [lines[k] for k in order]
        for a, b, v in lines:
            out.append("    C%d C%d %s" % (a, b, num(v)))
    out.append("ENDATA")
    return "\n".join(out) + "\n"


# ------------------------------------------------------------------------------------------------ building blocks
def dy(rs, size, k=24):
    """dyadic values in multiples of 1/8 (exact products / sums in double precision)."""
    return rs.randint(-k, k + 1, size) / 8.0


def int_psd(rs, n, rank=None, diag=0, kmax=3):
    """Integer PSD matrix B'B + diag*I (exactly symmetric, exactly PSD)."""
    B = rs.randint(-kmax, kmax + 1, (rank if rank is not None else n, n)).astype(float)
    return B.T @ B + diag * np.eye(n)


def sparse_dy(rs, m, n, dens, k=24):
    A = dy(rs, (m, n), k) * (rs.random_sample((m, n)) < dens)
    return A


def rows_through(rs, A, x0, kinds=("E", "L", "G", "R"), slack_zero=0.4):
    """Row bounds from x0 so that x0 is feasible; slack zero (active) with probability slack_zero."""
    m = A.shape[0]
    act = A @ x0
    rlo, rup = np.full(m, -INF), np.full(m, INF)
    for i in range(m):
        k = kinds[rs.randint(len(kinds))]
        s1 = 0.0 if rs.random_sample() < slack_zero else rs.randint(1, 17) / 8.0
        s2 = 0.0 if rs.random_sample() < slack_zero else rs.randint(1, 17) / 8.0
        if k == "E":
            rlo[i] = rup[i] = act[i]
        elif k == "L":
            rup[i] = act[i] + s1
        elif k == "G":
            rlo[i] = act[i] - s1
        else:
            rlo[i], rup[i] = act[i] - s1, act[i] + s2
    return rlo, rup


def box_through(rs, x0, kinds=("box", "box", "lo", "up", "free", "fixed"), slack_zero=0.35):
    n = len(x0)
    lo, up = np.full(n, -INF), np.full(n, INF)
    for j in range(n):
        k = kinds[rs.randint(len(kinds))]
        s1 = 0.0 if rs.random_sample() < slack_zero else rs.randint(1, 25) / 8.0
        s2 = 0.0 if rs.random_sample() < slack_zero else rs.randint(1, 25) / 8.0
        if k == "box":
            lo[j], up[j] = x0[j] - s1, x0[j] + s2
        elif k == "lo":
            lo[j] = x0[j] - s1
        elif k == "up":
            up[j] = x0[j] + s2
        elif k == "fixed":
            lo[j] = up[j] = x0[j]
        # free: both infinite
    return lo, up


def kkt_cost(rs, Qm, A, rlo, rup, lo, up, x0, zero_frac=0.3, ymax=5):
    """Cost c (minimisation sense) making x0 an optimal point of  min c'x + 0.5 x'Qm x  over the polyhedron.

    The multipliers are drawn with the right signs on the constraints active at x0 (exact test: act == bound);
    zero_frac of them are set to 0 (weakly active). KKT is sufficient for a convex Qm, so x0 is optimal.
    """
    m, n = A.shape
    act = A @ x0
    y = np.zeros(m)
    for i in range(m):
        al, au = rlo[i] > -INF and act[i] == rlo[i], rup[i] < INF and act[i] == rup[i]
        if not (al or au) or rs.random_sample() < zero_frac:
            continue
        mag = rs.randint(1, 8 * ymax) / 8.0
        y[i] = mag * (rs.choice([-1, 1]) if (al and au) else (1 if al else -1))
    z = np.zeros(n)
    for j in range(n):
        al, au = lo[j] > -INF and x0[j] == lo[j], up[j] < INF and x0[j] == up[j]
        if not (al or au) or rs.random_sample() < zero_frac:
            continue
        mag = rs.randint(1, 8 * ymax) / 8.0
        z[j] = mag * (rs.choice([-1, 1]) if (al and au) else (1 if al else -1))
    return A.T @ y + z - Qm @ x0


def finish_constructed(case, x0):
    case.x0 = np.asarray(x0, float)
    case.expect_status = "optimal"
    case.expect_obj = case.obj(case.x0)
    return case


def maybe_max(rs, Qm, c, p=0.3):
    """Present a minimisation problem as a max one (negate Q and c) with probability p."""
    if rs.random_sample() < p:
        return -Qm, -c, True
    return Qm, c, False


def constructed(cat, seed, rs, n, m, Q, dens=0.6, bk=None, rk=None, zero_frac=0.3, max_p=0.3, k=24, **kw):
    """Random polyhedron through a dyadic x0, KKT cost, random max presentation. Returns a finished QCase."""
    x0 = dy(rs, n, 16)
    A = sparse_dy(rs, m, n, dens, k) if m else np.zeros((0, n))
    rlo, rup = rows_through(rs, A, x0, **({"kinds": rk} if rk else {})) if m else (np.zeros(0), np.zeros(0))
    lo, up = box_through(rs, x0, **({"kinds": bk} if bk else {}))
    c = kkt_cost(rs, Q, A, rlo, rup, lo, up, x0, zero_frac)
    Qm, cm, mx = maybe_max(rs, Q, c, max_p)
    case = QCase(cat, seed, Qm, cm, A, rlo, rup, lo, up, maximize=mx, rowstyle=int(rs.randint(3)), **kw)
    return finish_constructed(case, x0)


def N(rs, lo, hi):
    return int(rs.randint(lo, hi + 1))


# ------------------------------------------------------------------------------------------------ categories
def g_pd_baseline(seed):
    rs = np.random.RandomState(seed)
    n = N(rs, 2, 18)
    return constructed("pd_baseline", seed, rs, n, N(rs, 0, n + 3), int_psd(rs, n, diag=1))


def g_psd_rank_def(seed):
    """PSD Q of low rank (including rank 0 and 1), all columns boxed so a minimiser exists."""
    rs = np.random.RandomState(seed)
    n = N(rs, 3, 20)
    Q = int_psd(rs, n, rank=N(rs, 0, max(1, n // 2)))
    bk = ("box", "box", "box", "fixed")
    return constructed("psd_rank_def", seed, rs, n, N(rs, 0, n + 2), Q, bk=bk)


def g_psd_free_null(seed):
    """PSD Q, free columns, equality rows only; c in range(Q)+range(A'), so an (unbounded-domain) minimiser exists."""
    rs = np.random.RandomState(seed)
    n = N(rs, 3, 14)
    Q = int_psd(rs, n, rank=N(rs, 1, n - 1))
    m = N(rs, 0, max(1, n // 2))
    x0 = dy(rs, n, 16)
    A = sparse_dy(rs, m, n, 0.7, 8) if m else np.zeros((0, n))
    b = A @ x0
    c = Q @ rs.randint(-3, 4, n).astype(float) + A.T @ rs.randint(-3, 4, m).astype(float)
    case = QCase("psd_free_null", seed, Q, c, A, b, b, np.full(n, -INF), np.full(n, INF), rowstyle=int(rs.randint(3)))
    case.expect_status = "optimal"  # c in range(Q)+range(A'): bounded below; the minimiser is not unique
    case.note = "free columns, equality rows, c in range(Q)+range(A')"
    return case


def g_unbounded(seed):
    """Feasible, with a verified recession direction d: Ad in the row cone, d in the column cone, Qd=0, c'd<0."""
    rs = np.random.RandomState(seed)
    kind = seed % 3  # 0 unit ray outside Q's support; 1 two-column ray; 2 ray inside the null space of a rank-def Q
    n = N(rs, 3, 12)
    p = N(rs, 1, n - 2)
    Q = np.zeros((n, n))
    Q[:p, :p] = int_psd(rs, p, diag=1)
    d = np.zeros(n)
    if kind == 0:
        d[p] = 1
    elif kind == 1:
        d[p], d[p + 1] = 1, 1
    else:
        # rank-deficient block: Q = B'B with B having a zero column pattern -> d = e_{p} is in null, add Q rows to it
        d[p] = 1
        Q[p + 1, p + 1] = 0
    x0 = dy(rs, n, 12)
    m = N(rs, 0, 4)
    A = np.round(rs.randint(-4, 5, (m, n)).astype(float))
    if kind == 1 and m:
        A[:, p + 1] = -A[:, p]  # A d = 0
    else:
        A[:, p] = 0
    lo, up = x0 - rs.randint(1, 8, n), x0 + rs.randint(1, 8, n)
    for j in range(n):
        if d[j] > 0:
            up[j] = INF
    rlo, rup = rows_through(rs, A, x0) if m else (np.zeros(0), np.zeros(0))
    c = rs.randint(-4, 5, n).astype(float)
    c[p] = -float(rs.randint(1, 6))
    if kind == 1:
        c[p + 1] = -float(rs.randint(1, 6))
    mx = rs.random_sample() < 0.3
    Qm, cm = (-Q, -c) if mx else (Q, c)
    case = QCase("unbounded", seed, Qm, cm, A, rlo, rup, lo, up, maximize=mx, rowstyle=int(rs.randint(3)))
    case.expect_status, case.ray, case.x0 = "unbounded", d, x0
    return case


def g_infeasible(seed):
    """Provably infeasible (integer data, margin >= 1e-3), PSD or PD Q."""
    rs = np.random.RandomState(seed)
    kind = seed % 5
    n = N(rs, 2, 12)
    Q = int_psd(rs, n, rank=N(rs, 0, n), diag=int(rs.randint(0, 2)))
    x0 = dy(rs, n, 8)
    A = np.round(rs.randint(-3, 4, (N(rs, 1, 4), n)).astype(float))
    rlo, rup = rows_through(rs, A, x0)
    lo, up = box_through(rs, x0, kinds=("box", "lo", "free"), slack_zero=0.2)
    a = rs.randint(1, 4, n).astype(float)
    gap = [1e-3, 1e-2, 0.5, 4.0][seed % 4]
    if kind == 0:  # a'x >= beyond the box maximum (all boxed)
        lo, up = x0 - 3, x0 + 3
        A = np.vstack([A, a]); rlo = np.append(rlo, a @ up + gap); rup = np.append(rup, INF)
    elif kind == 1:  # parallel contradictory rows
        A = np.vstack([A, a, a]); rlo = np.append(rlo, [a @ x0 + gap, -INF]); rup = np.append(rup, [INF, a @ x0])
    elif kind == 2:  # inconsistent dependent equalities: row3 = row1 + row2 with the wrong right-hand side
        r1, r2 = np.round(rs.randint(-3, 4, n)).astype(float), np.round(rs.randint(-3, 4, n)).astype(float)
        b1, b2 = float(r1 @ x0), float(r2 @ x0)
        A = np.vstack([A, r1, r2, r1 + r2]); rlo = np.append(rlo, [b1, b2, b1 + b2 + gap]); rup = np.append(rup, [b1, b2, b1 + b2 + gap])
    elif kind == 3:  # crossed column bounds on one column
        j = int(rs.randint(n)); lo = lo.copy(); up = up.copy(); lo[j], up[j] = x0[j] + gap, x0[j]
    else:  # all columns fixed at x0 but an equality row asks for a'x = a'x0 + gap
        lo, up = x0.copy(), x0.copy()
        A = np.vstack([A, a]); rlo = np.append(rlo, a @ x0 + gap); rup = np.append(rup, a @ x0 + gap)
    mx = rs.random_sample() < 0.3
    Qm, c = (-Q, -rs.randint(-3, 4, n).astype(float)) if mx else (Q, rs.randint(-3, 4, n).astype(float))
    case = QCase("infeasible", seed, Qm, c, A, rlo, rup, lo, up, maximize=mx, rowstyle=int(rs.randint(3)))
    case.expect_status = "infeasible"
    case.note = "kind %d gap %g" % (kind, gap)
    return case


def g_indefinite(seed):
    """Clearly indefinite or concave Q (min sense), bounded domain: a global-optimality claim is not defensible."""
    rs = np.random.RandomState(seed)
    n = N(rs, 2, 12)
    kind = seed % 4
    Qm = int_psd(rs, n, diag=1)
    if kind == 0:  # one flipped eigen-direction
        w, U = np.linalg.eigh(Qm); w[0] = -abs(w[0]) - 1.0; Qm = (U * w) @ U.T; Qm = np.round(Qm * 8) / 8
    elif kind == 1:  # negative definite
        Qm = -Qm
    elif kind == 2:  # zero diagonal with off-diagonal coupling: indefinite saddle
        Qm = np.zeros((n, n))
        for i in range(n):
            for j in range(i):
                Qm[i, j] = Qm[j, i] = float(rs.randint(-2, 3))
        Qm[0, 1] = Qm[1, 0] = 2.0
    else:  # PSD plus a small negative direction (-1e-3 .. -1)
        v = rs.randint(-2, 3, n).astype(float); v[0] = 1.0
        Qm = int_psd(rs, n, rank=max(1, n // 2)) - [1e-3, 1e-2, 0.1, 1.0][int(rs.randint(4))] * np.outer(v, v) / max(1.0, v @ v)
    Qm = (Qm + Qm.T) / 2
    x0 = dy(rs, n, 8)
    A = sparse_dy(rs, N(rs, 0, 3), n, 0.6, 8)
    rlo, rup = rows_through(rs, A, x0) if A.shape[0] else (np.zeros(0), np.zeros(0))
    lo, up = x0 - rs.randint(1, 6, n), x0 + rs.randint(1, 6, n)
    c = rs.randint(-4, 5, n).astype(float)
    mx = rs.random_sample() < 0.3
    Q, c = (-Qm, -c) if mx else (Qm, c)
    case = QCase("indefinite", seed, Q, c, A, rlo, rup, lo, up, maximize=mx, rowstyle=int(rs.randint(3)))
    case.x0 = x0
    case.note = "kind %d" % kind
    return case


def g_indef_convex_feasible(seed):
    """Q indefinite but PD on the null space of the equality rows: the problem is convex on its feasible set.
    Free columns, equality rows only. Reference: the KKT linear system (numpy), no solver."""
    rs = np.random.RandomState(seed)
    n = N(rs, 3, 10)
    m = N(rs, 1, n - 1)
    A = np.round(rs.randint(-3, 4, (m, n)).astype(float))
    while np.linalg.matrix_rank(A) < m:
        A = np.round(rs.randint(-3, 4, (m, n)).astype(float))
    U, s, Vt = np.linalg.svd(A)
    Z = Vt[m:].T  # null-space basis
    x0 = dy(rs, n, 12)
    # indefinite Q with Z'QZ PD: Q = Z P Z' + (A' S + S' A)-type terms with big negative curvature along row space
    P = int_psd(rs, n - m, diag=1)
    R = Vt[:m].T
    Q = Z @ P @ Z.T - 3.0 * (R @ R.T) + np.round(rs.randn(n, n) * 0.0)
    Q = np.round((Q + Q.T) / 2 * 64) / 64  # dyadic, still indefinite (A-rows directions negative)
    if np.linalg.eigvalsh(Z.T @ Q @ Z)[0] <= 1e-3 or np.linalg.eigvalsh(Q)[0] >= -1e-3:
        Q = Z @ P @ Z.T - 3.0 * (R @ R.T)
        Q = (Q + Q.T) / 2
    c = rs.randint(-4, 5, n).astype(float)
    b = A @ x0
    case = QCase("indef_convex_feasible", seed, Q, c, A, b, b, np.full(n, -INF), np.full(n, INF), rowstyle=int(rs.randint(3)))
    case.x0 = x0
    case.expect_status = "optimal"
    # exact-ish optimum from the KKT system  [Q A'; A 0][x; -y] = [-c; b]
    K = np.block([[Q, A.T], [A, np.zeros((m, m))]])
    sol = np.linalg.solve(K, np.concatenate([-c, b]))
    case.x0, case.expect_obj = sol[:n], case.obj(sol[:n])
    case.note = "convex on {Ax=b}, Q indefinite, lam_min(Q)=%.3g" % np.linalg.eigvalsh(Q)[0]
    return case


def g_near_psd(seed):
    """PSD Q plus eps * vv' with eps from -1e-4 to +1e-14 (relative to |Q|): the convexity boundary."""
    rs = np.random.RandomState(seed)
    n = N(rs, 2, 12)
    Q = int_psd(rs, n, rank=N(rs, 1, n))
    v = rs.randn(n)
    # put the perturbation in a null direction of Q when there is one, else along an arbitrary direction
    w, U = np.linalg.eigh(Q)
    v = U[:, 0] if w[0] < 1e-9 else v / np.linalg.norm(v)
    eps = [-1e-4, -1e-6, -1e-8, -1e-10, -1e-12, -1e-14, 1e-14, 1e-10][seed % 8] * max(1.0, w[-1])
    Qm = (Q + eps * np.outer(v, v))
    Qm = (Qm + Qm.T) / 2
    x0 = dy(rs, n, 8)
    A = sparse_dy(rs, N(rs, 0, 3), n, 0.6, 8)
    rlo, rup = rows_through(rs, A, x0) if A.shape[0] else (np.zeros(0), np.zeros(0))
    lo, up = x0 - rs.randint(1, 6, n), x0 + rs.randint(1, 6, n)
    c = rs.randint(-4, 5, n).astype(float)
    case = QCase("near_psd", seed, Qm, c, A, rlo, rup, lo, up, rowstyle=int(rs.randint(3)))
    case.x0 = x0
    case.note = "eps=%g" % eps
    return case


def g_zero_q(seed):
    """Q = 0 written as an empty QUADOBJ, explicit zeros, cancelling pairs, or no section: an LP through the QP path."""
    rs = np.random.RandomState(seed)
    n = N(rs, 2, 14)
    case = constructed("zero_q", seed, rs, n, N(rs, 1, n + 2), np.zeros((n, n)), zero_frac=0.2)
    kind = seed % 4
    if kind == 0:
        case.qempty = True
    elif kind == 1:
        case.qzeros = True
    elif kind == 2:  # tiny nonzero Q entry on a fixed column only: effectively zero Hessian
        j = int(np.flatnonzero(case.lo == case.up)[0]) if np.any(case.lo == case.up) else 0
        case.Q = case.Q.copy(); case.Q[j, j] = 1e-12 * (-1 if case.maximize else 1)
        case.lo[j] = case.up[j] = case.x0[j]
        case.expect_obj = case.obj(case.x0)
        case.note = "1e-12 curvature on column %d" % j
    return case


def g_degenerate_cons(seed):
    """Duplicate / dependent / parallel / zero rows and a degenerate vertex (more active constraints than columns)."""
    rs = np.random.RandomState(seed)
    n = N(rs, 2, 10)
    Q = int_psd(rs, n, rank=N(rs, 0, n), diag=int(rs.randint(0, 2)))
    x0 = dy(rs, n, 8)
    base = sparse_dy(rs, N(rs, 1, n), n, 0.7, 6)
    rows = [r for r in base]
    kind = seed % 5
    for _ in range(N(rs, 1, 4)):
        r = rows[rs.randint(len(rows))]
        if kind == 0:
            rows.append(r.copy())                               # exact duplicate
        elif kind == 1:
            rows.append(r * float(rs.choice([-2, 2, 0.5, 8])))  # scaled duplicate
        elif kind == 2:
            rows.append(r + rows[rs.randint(len(rows))])        # dependent (sum of two)
        elif kind == 3:
            rows.append(np.zeros(n))                            # zero row
        else:
            rows.append(r.copy())
    extra = n + 2 - len(rows)  # degenerate vertex: pad with rows all active at x0
    for _ in range(max(0, extra)):
        rows.append(sparse_dy(rs, 1, n, 0.7, 6)[0])
    A = np.array(rows)
    act = A @ x0
    m = A.shape[0]
    rlo, rup = np.full(m, -INF), np.full(m, INF)
    for i in range(m):
        t = rs.randint(4)
        if t == 0:
            rlo[i] = rup[i] = act[i]
        elif t == 1:
            rup[i] = act[i]
        elif t == 2:
            rlo[i] = act[i]
        else:
            rlo[i], rup[i] = act[i] - 1, act[i] + 1
    if kind == 4:  # parallel rows with a looser copy
        r = rows[0]; A = np.vstack([A, r]); rlo = np.append(rlo, -INF); rup = np.append(rup, float(r @ x0) + 1)
    lo, up = box_through(rs, x0, kinds=("box", "lo", "up", "free"), slack_zero=0.6)
    c = kkt_cost(rs, Q, A, rlo, rup, lo, up, x0, zero_frac=0.4)
    Qm, cm, mx = maybe_max(rs, Q, c)
    case = QCase("degenerate_cons", seed, Qm, cm, A, rlo, rup, lo, up, maximize=mx, rowstyle=int(rs.randint(3)))
    case.note = "kind %d" % kind
    return finish_constructed(case, x0)


def g_active_set(seed):
    """Optimum on many bounds and rows, strictly and weakly active multipliers mixed."""
    rs = np.random.RandomState(seed)
    n = N(rs, 4, 30)
    Q = int_psd(rs, n, rank=N(rs, 1, max(1, n // 3)), diag=int(rs.randint(0, 2)))
    return constructed("active_set", seed, rs, n, N(rs, 2, n), Q, zero_frac=float(rs.choice([0.0, 0.3, 0.8])),
                       bk=("box", "box", "fixed", "lo", "up"), rk=("E", "L", "G", "R"))


def g_bounds_special(seed):
    """Fixed, free, negative lower, upper-only (MI) and zero-width ranges; Q couples the special columns."""
    rs = np.random.RandomState(seed)
    n = N(rs, 3, 14)
    Q = int_psd(rs, n, rank=N(rs, 1, n), diag=int(rs.randint(0, 2)))
    kind = seed % 3
    bk = [("fixed", "fixed", "box", "lo"), ("free", "free", "box", "up"), ("box", "lo", "up", "fixed", "free")][kind]
    case = constructed("bounds_special", seed, rs, n, N(rs, 0, n), Q, bk=bk)
    # negative lower bounds and negative upper-only bounds through a shifted x0 are already produced by dy()
    return case


def g_bounds_huge(seed):
    """Inactive bounds of size 1e6..1e19 (finite) and 1e20..1e30 (infinite by convention) on columns and rows."""
    rs = np.random.RandomState(seed)
    n = N(rs, 2, 12)
    Q = int_psd(rs, n, rank=N(rs, 1, n), diag=int(rs.randint(0, 2)))
    case = constructed("bounds_huge", seed, rs, n, N(rs, 0, n), Q, max_p=0.2)
    mag = [1e6, 1e9, 1e12, 1e15, 1e19, 1e20, 1e25, 1e30][seed % 8]
    x0 = case.x0
    for j in range(n):  # move inactive (strictly slack) bounds out to +-mag
        if case.lo[j] > -INF and x0[j] > case.lo[j] and rs.random_sample() < 0.6:
            case.lo[j] = -mag
        if case.up[j] < INF and x0[j] < case.up[j] and rs.random_sample() < 0.6:
            case.up[j] = mag
    act = case.A @ x0 if case.m else np.zeros(0)
    for i in range(case.m):
        if case.rlo[i] > -INF and act[i] > case.rlo[i] and rs.random_sample() < 0.6:
            case.rlo[i] = -mag
        if case.rup[i] < INF and act[i] < case.rup[i] and rs.random_sample() < 0.6:
            case.rup[i] = mag
    for a in ("rlo", "rup", "lo", "up"):
        v = getattr(case, a)
        setattr(case, a, np.where(v >= 1e20, INF, np.where(v <= -1e20, -INF, v)))
    case.note = "magnitude %g" % mag
    return case


def g_bounds_tight(seed):
    """Intervals of width 1e-9 .. 1e-3 around x0 on columns and rows (almost fixed)."""
    rs = np.random.RandomState(seed)
    n = N(rs, 2, 12)
    Q = int_psd(rs, n, rank=N(rs, 1, n), diag=int(rs.randint(0, 2)))
    case = constructed("bounds_tight", seed, rs, n, N(rs, 0, n), Q)
    w = [1e-9, 1e-7, 1e-5, 1e-3][seed % 4]
    sel = rs.random_sample(n) < 0.5
    case.lo = np.where(sel & (case.lo > -INF) & (case.lo < case.x0), case.x0 - w, case.lo)
    case.up = np.where(sel & (case.up < INF) & (case.up > case.x0), case.x0 + w, case.up)
    case.note = "width %g" % w
    # x0 stays feasible and optimal: only inactive bounds moved closer (active ones are unchanged)
    return case


def g_scale_q(seed):
    """Hessian scaled by 1e-6..1e6 against fixed constraints (cost rebuilt so the optimum is known)."""
    rs = np.random.RandomState(seed)
    n = N(rs, 2, 14)
    s = 10.0 ** [-6, -4, -2, 2, 4, 6][seed % 6]
    Q = s * int_psd(rs, n, rank=N(rs, 1, n), diag=int(rs.randint(0, 2)))
    case = constructed("scale_q", seed, rs, n, N(rs, 0, n), Q)
    case.note = "Q scale %g" % s
    return case


def g_scale_cols(seed):
    """x = D x' with D = 10^U(-4,4): columns of A, Q, c, the bounds and x0 are rescaled together."""
    rs = np.random.RandomState(seed)
    n = N(rs, 2, 14)
    Q = int_psd(rs, n, rank=N(rs, 1, n), diag=1)
    b = constructed("scale_cols", seed, rs, n, N(rs, 0, n), Q, max_p=0.0)
    D = 10.0 ** rs.uniform(-4, 4, n)
    case = QCase("scale_cols", seed, (b.Q * D[:, None]) * D[None, :], b.c * D, b.A * D[None, :], b.rlo, b.rup,
                 b.lo / D, b.up / D, rowstyle=b.rowstyle)
    case.note = "column scale span %.1e" % (D.max() / D.min())
    return finish_constructed(case, b.x0 / D)


def g_scale_rows(seed):
    """Rows multiplied by 10^U(-6,6) (bounds with them)."""
    rs = np.random.RandomState(seed)
    n = N(rs, 2, 14)
    Q = int_psd(rs, n, rank=N(rs, 1, n), diag=1)
    b = constructed("scale_rows", seed, rs, n, N(rs, 1, n + 1), Q, max_p=0.0)
    S = 10.0 ** rs.uniform(-6, 6, b.m)
    # the dual scales inversely, cost stays valid: A'y with y' = y/S is the same cost vector
    case = QCase("scale_rows", seed, b.Q, b.c, b.A * S[:, None], b.rlo * S, b.rup * S, b.lo, b.up, rowstyle=b.rowstyle)
    case.note = "row scale span %.1e" % (S.max() / S.min())
    return finish_constructed(case, b.x0)


def g_ill_cond(seed):
    """Q with condition number 1e6 .. 1e12 (eigenvalues log-uniform); box-bounded; reference-only (no exact optimum)."""
    rs = np.random.RandomState(seed)
    n = N(rs, 3, 14)
    U, _ = np.linalg.qr(rs.randn(n, n))
    kappa = 10.0 ** [6, 8, 10, 12][seed % 4]
    ev = np.exp(rs.uniform(0, math.log(kappa), n)); ev[0] = 1.0; ev[-1] = kappa
    Q = (U * ev) @ U.T
    Q = (Q + Q.T) / 2
    x0 = dy(rs, n, 8)
    A = sparse_dy(rs, N(rs, 0, 4), n, 0.6, 8)
    rlo, rup = rows_through(rs, A, x0) if A.shape[0] else (np.zeros(0), np.zeros(0))
    lo, up = x0 - rs.randint(1, 8, n), x0 + rs.randint(1, 8, n)
    c = rs.randn(n) * 10
    case = QCase("ill_cond", seed, Q, c, A, rlo, rup, lo, up, rowstyle=int(rs.randint(3)))
    case.x0 = x0
    case.expect_status = "optimal"
    case.note = "kappa %g" % kappa
    return case


def g_maximize(seed):
    """OBJSENSE MAX with concave Q (convex), constants of both signs; half also with a free-ish row mix."""
    rs = np.random.RandomState(seed)
    n = N(rs, 2, 14)
    Qm = int_psd(rs, n, rank=N(rs, 1, n), diag=int(rs.randint(0, 2)))
    case = constructed("maximize", seed, rs, n, N(rs, 0, n + 1), Qm, max_p=1.0)
    case.const = float(rs.choice([-1e4, -3.5, 0.0, 2.25, 1e6]))
    case.expect_obj = case.obj(case.x0)
    return case


def g_obj_const(seed):
    rs = np.random.RandomState(seed)
    n = N(rs, 2, 10)
    case = constructed("obj_const", seed, rs, n, N(rs, 0, n), int_psd(rs, n, rank=N(rs, 0, n), diag=int(rs.randint(0, 2))),
                       max_p=0.5)
    case.const = float(rs.choice([-1e9, -1e3, -0.125, 0.125, 1e3, 1e9, 123456.789]))
    case.expect_obj = case.obj(case.x0)
    return case


def g_qform(seed):
    """The same Q written as one triangle (either token order), QMATRIX, QSECTION, split duplicates, zero filler."""
    rs = np.random.RandomState(seed)
    n = N(rs, 2, 10)
    case = constructed("qform", seed, rs, n, N(rs, 0, n), int_psd(rs, n, rank=N(rs, 1, n), diag=int(rs.randint(0, 2))),
                       max_p=0.3)
    kind = seed % 6
    case.qform = ["lower", "upper", "full", "qsection", "lower", "lower"][kind]
    case.qdups = kind == 4
    case.qzeros = kind == 5
    case.note = ["QUADOBJ lower", "QUADOBJ swapped tokens", "QMATRIX", "QSECTION OBJ", "split duplicates", "zero filler"][kind]
    return case


def g_empty_structure(seed):
    rs = np.random.RandomState(seed)
    kind = seed % 7
    n = N(rs, 1, 8)
    Q = int_psd(rs, n, rank=N(rs, 1, n), diag=1)
    if kind == 0:  # no rows
        return constructed("empty_structure", seed, rs, n, 0, Q)
    if kind == 1:  # a column in no row and with no cost/curvature
        Q[-1, :] = 0; Q[:, -1] = 0
        n_ = n
        case = constructed("empty_structure", seed, rs, n_, N(rs, 0, n_), Q)
        if case.m:
            case.A[:, -1] = 0
        case.c[-1] = 0
        case.lo[-1], case.up[-1] = -1.0, 1.0
        case.x0[-1] = 0.0
        case.rlo, case.rup = rows_through(rs, case.A, case.x0) if case.m else (case.rlo, case.rup)
        case.c = kkt_cost(rs, case.qmin(), case.A, case.rlo, case.rup, case.lo, case.up, case.x0) * (-1 if case.maximize else 1)
        case.expect_obj = case.obj(case.x0)
        return case
    if kind == 2:  # n = 1
        return constructed("empty_structure", seed, rs, 1, N(rs, 0, 2), np.array([[float(rs.randint(1, 5))]]))
    if kind == 3:  # all columns fixed
        return constructed("empty_structure", seed, rs, n, N(rs, 0, 4), Q, bk=("fixed",))
    if kind == 4:  # row with no coefficients, feasible (0 inside its range)
        case = constructed("empty_structure", seed, rs, n, N(rs, 1, 3), Q)
        case.A = np.vstack([case.A, np.zeros(case.n)])
        case.rlo = np.append(case.rlo, -1.0); case.rup = np.append(case.rup, 2.0)
        return case
    if kind == 5:  # row with no coefficients, infeasible (0 outside its range)
        case = constructed("empty_structure", seed, rs, n, N(rs, 0, 3), Q)
        case.A = np.vstack([case.A, np.zeros(case.n)]) if case.m else np.zeros((1, case.n))
        case.rlo = np.append(case.rlo, 1e-3); case.rup = np.append(case.rup, INF)
        case.expect_status, case.expect_obj, case.x0 = "infeasible", None, None
        return case
    # kind 6: column present only in Q (no row entries), others in rows
    case = constructed("empty_structure", seed, rs, n, N(rs, 1, 4), Q)
    if case.n > 1:
        case.A[:, 0] = 0
        act = case.A @ case.x0
        case.rlo, case.rup = rows_through(rs, case.A, case.x0)
        case.c = kkt_cost(rs, case.qmin(), case.A, case.rlo, case.rup, case.lo, case.up, case.x0) * (-1 if case.maximize else 1)
        case.expect_obj = case.obj(case.x0)
    return case


def g_portfolio(seed):
    """Factor-model Markowitz: x >= 0, sum x = 1, return floor; Q = F'F + diag; reference-only."""
    rs = np.random.RandomState(seed)
    n = N(rs, 5, 60)
    k = N(rs, 1, 5)
    F = rs.randn(k, n) * 0.2
    Q = 2 * (F.T @ F + np.diag(rs.uniform(0.01, 0.2, n)))
    Q = (Q + Q.T) / 2
    mu = rs.uniform(0.0, 0.2, n)
    A = np.vstack([np.ones(n), mu])
    rlo = np.array([1.0, float(np.quantile(mu, 0.5))])
    rup = np.array([1.0, INF])
    case = QCase("portfolio", seed, Q, np.zeros(n), A, rlo, rup, np.zeros(n), np.full(n, INF if seed % 2 else 0.5),
                 rowstyle=int(rs.randint(3)))
    case.expect_status = "optimal"
    return case


def g_large_sparse(seed):
    """n = 150..400, banded/sparse PSD Q, sparse rows; KKT-constructed optimum."""
    rs = np.random.RandomState(seed)
    n = N(rs, 150, 400)
    B = (rs.randint(-3, 4, (n // 2, n)) * (rs.random_sample((n // 2, n)) < 0.02)).astype(float)
    Q = B.T @ B + np.diag((rs.random_sample(n) < 0.5).astype(float))
    return constructed("large_sparse", seed, rs, n, n // 2, Q, dens=0.03, max_p=0.2)


def g_fuzz_mixed(seed):
    """Random mix: rank, bounds, rows, scaling, form and sense chosen at random; truth by KKT construction."""
    rs = np.random.RandomState(seed)
    n = N(rs, 2, 16)
    Q = int_psd(rs, n, rank=N(rs, 0, n), diag=int(rs.randint(0, 2)))
    case = constructed("fuzz_mixed", seed, rs, n, N(rs, 0, n + 3), Q, dens=float(rs.choice([0.2, 0.6, 1.0])),
                       zero_frac=float(rs.choice([0.0, 0.3, 0.7])), max_p=0.4)
    case.qform = str(rs.choice(["lower", "lower", "upper", "full"]))
    if rs.random_sample() < 0.3:
        case.const = float(rs.randint(-50, 50))
        case.expect_obj = case.obj(case.x0)
    return case


# (name, generator, count) in category order; kind k in the seed base is the index here
ALL = [
    ("pd_baseline", g_pd_baseline, 80),
    ("psd_rank_def", g_psd_rank_def, 80),
    ("psd_free_null", g_psd_free_null, 60),
    ("unbounded", g_unbounded, 90),
    ("infeasible", g_infeasible, 100),
    ("indefinite", g_indefinite, 80),
    ("indef_convex_feasible", g_indef_convex_feasible, 40),
    ("near_psd", g_near_psd, 80),
    ("zero_q", g_zero_q, 60),
    ("degenerate_cons", g_degenerate_cons, 100),
    ("active_set", g_active_set, 80),
    ("bounds_special", g_bounds_special, 90),
    ("bounds_huge", g_bounds_huge, 80),
    ("bounds_tight", g_bounds_tight, 60),
    ("scale_q", g_scale_q, 60),
    ("scale_cols", g_scale_cols, 60),
    ("scale_rows", g_scale_rows, 60),
    ("ill_cond", g_ill_cond, 60),
    ("maximize", g_maximize, 60),
    ("obj_const", g_obj_const, 60),
    ("qform", g_qform, 60),
    ("empty_structure", g_empty_structure, 70),
    ("portfolio", g_portfolio, 40),
    ("large_sparse", g_large_sparse, 20),
    ("fuzz_mixed", g_fuzz_mixed, 150),
]
BASE = {name: SEED_BASE + 1000 * k for k, (name, _, _) in enumerate(ALL)}
GEN = {name: g for name, g, _ in ALL}
COUNT = {name: c for name, _, c in ALL}


def gen(cat, i):
    case = GEN[cat](BASE[cat] + i)
    assert case.cat == cat
    return case
