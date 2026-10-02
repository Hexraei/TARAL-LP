#!/usr/bin/env python3
"""Seeded randomized convex-QP suites for the C++ engine (QUADOBJ -> interior point), checked against HiGHS.

Objective is 0.5 x'Qx + c'x (MPS QUADOBJ, lower triangle). Each case is written to MPS, solved by the CLI and by
HiGHS (readModel on the same file); every reported solution is re-checked here against the generated data. Suites:
  a random      convex Q = B'B + diag, equality + coupling inequality + ranged rows, mixed bounds
  b rankdef     rank-deficient PSD Q (B wide)
  c illcond     eigenvalues spanning 1e-6..1e3
  d active      optimum on many bounds / rows (strongly pulling linear term)
  e infeasible  contradictory rows/bounds
  f unbounded   PSD Q with a free descent direction in its null space
  g nonconvex   indefinite Q: engine should say nonconvex; claiming optimal is reported (never "wrong" unless
                the point is infeasible, mis-reported, or worse than HiGHS' local result)
  h maximize    OBJSENSE MAX with negative-definite/semidefinite Q

Usage: ~/.venvs/taral-gpu/bin/python benchmarks/qp_tests.py [--engine out/taral] [--per-suite 25] [--seed 1]
Exit code 1 on any WRONG case.
"""
import argparse, json, math, os, subprocess, sys, tempfile
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import highspy

INF = math.inf


class Case:
    def __init__(self, name, Q, c, A, rlo, rup, lo, up, maximize=False, truth=None, nonconvex=False):
        self.name, self.Q, self.c, self.A = name, np.asarray(Q, float), np.asarray(c, float), np.asarray(A, float)
        self.rlo, self.rup = np.asarray(rlo, float), np.asarray(rup, float)
        self.lo, self.up = np.asarray(lo, float), np.asarray(up, float)
        self.maximize, self.truth, self.nonconvex = maximize, truth, nonconvex

    def obj(self, x):
        return float(self.c @ x + 0.5 * x @ self.Q @ x)


def num(v):
    return repr(float(v)) if v != int(v) else str(int(v))


def write_mps(case, path):
    m, n = case.A.shape
    out = ["NAME " + case.name]
    if case.maximize:
        out += ["OBJSENSE", "    MAX"]
    out += ["ROWS", " N obj"]
    kinds, rhs, rng = [], [], {}
    for i in range(m):
        lo, up = case.rlo[i], case.rup[i]
        if lo == up:
            kinds.append("E"); rhs.append(lo)
        elif lo == -INF:
            kinds.append("L"); rhs.append(up)
        elif up == INF:
            kinds.append("G"); rhs.append(lo)
        else:
            kinds.append("L"); rhs.append(up); rng[i] = up - lo
        out.append(" %s r%d" % (kinds[-1], i))
    out.append("COLUMNS")
    for j in range(n):
        out.append("    c%d obj %s" % (j, num(case.c[j])))
        for i in range(m):
            if case.A[i, j] != 0:
                out.append("    c%d r%d %s" % (j, i, num(case.A[i, j])))
    out.append("RHS")
    for i in range(m):
        if rhs[i] != 0:
            out.append("    rhs r%d %s" % (i, num(rhs[i])))
    if rng:
        out.append("RANGES")
        for i, r in rng.items():
            out.append("    rng r%d %s" % (i, num(r)))
    out.append("BOUNDS")
    for j in range(n):
        lo, up = case.lo[j], case.up[j]
        if lo == up:
            out.append(" FX bnd c%d %s" % (j, num(lo)))
        elif lo == -INF and up == INF:
            out.append(" FR bnd c%d" % j)
        else:
            out.append((" MI bnd c%d" % j) if lo == -INF else (" LO bnd c%d %s" % (j, num(lo))))
            out.append((" PL bnd c%d" % j) if up == INF else (" UP bnd c%d %s" % (j, num(up))))
    out.append("QUADOBJ")
    for j in range(n):
        for i in range(j, n):
            if case.Q[i, j] != 0:
                out.append("    c%d c%d %s" % (j, i, num(case.Q[i, j])))
    out.append("ENDATA")
    with open(path, "w") as f:
        f.write("\n".join(out) + "\n")


def violation(case, x):
    act = case.A @ x
    w = 0.0
    for v, lo, up in list(zip(act, case.rlo, case.rup)) + list(zip(x, case.lo, case.up)):
        if lo > -INF:
            w = max(w, (lo - v) / (1 + abs(lo)))
        if up < INF:
            w = max(w, (v - up) / (1 + abs(up)))
    return w


def run_highs(path):
    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    h.setOptionValue("time_limit", 60.0)
    h.setOptionValue("threads", 1)
    h.readModel(path)
    h.run()
    return h.modelStatusToString(h.getModelStatus()), h.getInfo().objective_function_value


def run_taral(engine, path):
    js, sol = path + ".json", path + ".sol"
    for p in (js, sol):
        if os.path.exists(p):
            os.remove(p)
    subprocess.run([engine, path, "--time-limit", "20", "--json", js, "--sol", sol], capture_output=True, timeout=300)
    with open(js) as f:
        r = json.load(f)
    x = None
    if os.path.exists(sol):
        vals = dict(line.split() for line in open(sol))
        x = np.array([float(vals["c%d" % j]) for j in range(len(vals))])
    return r, x


def close(a, b, tol=1e-6):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def judge(case, r, x, hst, hobj):
    st = r["status"]
    sense = -1 if case.maximize else 1
    truth = case.truth or {"Optimal": "optimal", "Infeasible": "infeasible", "Unbounded": "unbounded"}.get(hst)
    if x is not None and r.get("has_solution", True) and st == "optimal":
        v = violation(case, x)
        if v > 1e-6:
            return "wrong", "reported point violates the model by %.3g" % v
        if not close(r["objective"], case.obj(x)):
            return "wrong", "objective %r does not match the point (%r)" % (r["objective"], case.obj(x))
    if case.nonconvex:
        if st == "nonconvex":
            return "ok", ""
        if st == "optimal" and x is not None:
            if hst == "Optimal" and sense * case.obj(x) > sense * hobj + 1e-6 * max(1, abs(hobj)):
                return "wrong", "optimal %r worse than HiGHS local %r" % (case.obj(x), hobj)
            return "undecided", "claimed optimal on indefinite Q (obj %r, HiGHS %s %r)" % (case.obj(x), hst, hobj)
        return "undecided", st
    if st == "optimal":
        if truth is None and hst not in ("Infeasible", "Unbounded", "Optimal"):
            return "undecided", "HiGHS gave no verdict (%s); point is feasible" % hst
        if truth != "optimal":
            return "wrong", "optimal but truth/HiGHS says %s" % (truth or hst)
        if x is not None and sense * case.obj(x) < sense * hobj - 1e-6 * max(1, abs(hobj)):
            return "undecided", "feasible point beats HiGHS optimum (%r vs %r): HiGHS suboptimal" % (case.obj(x), hobj)
        if x is None or not close(r["objective"], hobj):
            return "wrong", "objective %r vs HiGHS %r" % (r["objective"], hobj)
        return "ok", ""
    if st in ("infeasible", "unbounded", "dual_infeasible"):
        want = "unbounded" if st == "dual_infeasible" else st
        if truth == want:
            return "ok", ""
        if want == "unbounded" and hst == "Primal infeasible or unbounded" and truth is None:
            return "undecided", "HiGHS infeasible-or-unbounded"
        if want == "infeasible" and truth is None and hst == "Primal infeasible or unbounded":
            return "undecided", "HiGHS infeasible-or-unbounded"
        return "wrong", "%s but truth/HiGHS says %s" % (st, truth or hst)
    if st == "nonconvex":
        return "wrong", "nonconvex reported for a convex model"
    return "undecided", st


# ---------------------------------------------------------------- generators
def psd(rng, n, rank=None, diag=0.0):
    B = rng.normal(size=(rank or n, n))
    return B.T @ B + diag * np.eye(n)


def rows(rng, m, n, x0, ranged=True):
    """Coupling rows through x0: mix of E, L, G, ranged."""
    A = rng.normal(size=(m, n)) * (rng.random((m, n)) < 0.6)
    act = A @ x0
    rlo, rup = np.full(m, -INF), np.full(m, INF)
    for i in range(m):
        t = rng.integers(0, 4 if ranged else 3)
        s = rng.uniform(0.1, 2)
        if t == 0: rlo[i] = rup[i] = act[i]
        elif t == 1: rup[i] = act[i] + s
        elif t == 2: rlo[i] = act[i] - s
        else: rlo[i], rup[i] = act[i] - s, act[i] + rng.uniform(0.1, 2)
    return A, rlo, rup


def boxes(rng, n, x0, free=0.15):
    lo, up = x0 - rng.uniform(0.2, 3, n), x0 + rng.uniform(0.2, 3, n)
    for j in range(n):
        r = rng.random()
        if r < free: lo[j], up[j] = -INF, INF
        elif r < 0.35: up[j] = INF
        elif r < 0.5: lo[j] = -INF
    return lo, up


def base(rng, k, name, Q, free=0.15, cscale=3.0):
    n = Q.shape[0]
    x0 = rng.normal(size=n)
    A, rlo, rup = rows(rng, int(rng.integers(2, n + 3)), n, x0)
    lo, up = boxes(rng, n, x0, free)
    c = rng.normal(size=n) * cscale
    return A, rlo, rup, lo, up, c


def g_random(rng, k):
    n = int(rng.integers(4, 25))
    Q = psd(rng, n, diag=0.05)
    A, rlo, rup, lo, up, c = base(rng, k, "", Q, free=0.0)
    return Case("random%d" % k, Q, c, A, rlo, rup, lo, up)


def g_rankdef(rng, k):
    n = int(rng.integers(6, 25))
    Q = psd(rng, n, rank=int(rng.integers(1, max(2, n // 2))))
    A, rlo, rup, lo, up, c = base(rng, k, "", Q, free=0.0)
    lo, up = np.where(np.isinf(lo), -5, lo), np.where(np.isinf(up), 5, up)  # keep bounded, no null-space escape
    return Case("rankdef%d" % k, Q, c, A, rlo, rup, lo, up)


def g_illcond(rng, k):
    n = int(rng.integers(4, 20))
    U, _ = np.linalg.qr(rng.normal(size=(n, n)))
    ev = 10 ** rng.uniform(-6, 3, n)
    ev[0], ev[-1] = 1e-6, 1e3
    Q = U @ np.diag(ev) @ U.T
    Q = (Q + Q.T) / 2
    A, rlo, rup, lo, up, c = base(rng, k, "", Q, free=0.0)
    lo, up = np.where(np.isinf(lo), -10, lo), np.where(np.isinf(up), 10, up)
    return Case("illcond%d" % k, Q, c, A, rlo, rup, lo, up)


def g_active(rng, k):
    n = int(rng.integers(8, 30))
    Q = psd(rng, n, rank=max(1, n // 3), diag=0.01)
    A, rlo, rup, lo, up, c = base(rng, k, "", Q, free=0.0, cscale=40.0)
    lo, up = np.where(np.isinf(lo), -4, lo), np.where(np.isinf(up), 4, up)
    return Case("active%d" % k, Q, c, A, rlo, rup, lo, up)


def g_infeasible(rng, k):
    n = int(rng.integers(3, 15))
    Q = psd(rng, n, diag=0.1)
    A, rlo, rup, lo, up, c = base(rng, k, "", Q, free=0.0)
    lo, up = np.where(np.isinf(lo), -5, lo), np.where(np.isinf(up), 5, up)
    a = np.abs(rng.normal(size=n)) + 0.1  # sum a_j x_j >= beyond the box maximum
    A = np.vstack([A, a])
    rlo = np.append(rlo, a @ up + rng.uniform(0.5, 3))
    rup = np.append(rup, INF)
    return Case("infeas%d" % k, Q, c, A, rlo, rup, lo, up, truth="infeasible")


def g_unbounded(rng, k):
    n = int(rng.integers(4, 15))
    p = int(rng.integers(1, n - 1))  # columns p.. are outside Q's support
    Q = np.zeros((n, n))
    Q[:p, :p] = psd(rng, p, diag=0.1)
    x0 = rng.normal(size=n)
    A, rlo, rup = rows(rng, int(rng.integers(0, 4)), n, x0, ranged=False)
    A[:, p] = 0  # direction e_p: leaves all rows unchanged
    lo, up = np.full(n, -3.0), np.full(n, 3.0)
    lo[p], up[p] = x0[p] - 1, INF
    c = rng.normal(size=n)
    c[p] = -rng.uniform(0.1, 5)
    return Case("unbd%d" % k, Q, c, A, rlo, rup, lo, up, truth="unbounded")


def g_nonconvex(rng, k):
    n = int(rng.integers(3, 15))
    Q = rng.normal(size=(n, n)); Q = (Q + Q.T) / 2
    Q -= (np.linalg.eigvalsh(Q)[0] - rng.uniform(0.5, 3)) * 0  # keep as is; force a negative eigenvalue below
    w, U = np.linalg.eigh(Q)
    w[0] = -abs(w[0]) - 1.0
    Q = (U * w) @ U.T
    Q = (Q + Q.T) / 2
    x0 = rng.normal(size=n)
    A, rlo, rup = rows(rng, int(rng.integers(1, n)), n, x0)
    lo, up = x0 - rng.uniform(0.5, 3, n), x0 + rng.uniform(0.5, 3, n)
    return Case("nonconvex%d" % k, Q, rng.normal(size=n), A, rlo, rup, lo, up, nonconvex=True)


def g_maximize(rng, k):
    n = int(rng.integers(4, 20))
    Qp = psd(rng, n, rank=int(rng.integers(max(1, n // 2), n + 1)), diag=0.05)
    A, rlo, rup, lo, up, c = base(rng, k, "", Qp, free=0.0)
    return Case("max%d" % k, -Qp, c, A, rlo, rup, lo, up, maximize=True)


SUITES = {"a_random": g_random, "b_rankdef": g_rankdef, "c_illcond": g_illcond, "d_active": g_active,
          "e_infeasible": g_infeasible, "f_unbounded": g_unbounded, "g_nonconvex": g_nonconvex, "h_maximize": g_maximize}


def tiny(engine):
    """min x + 0.5*2*x^2, x >= 1 -> 2 (checks the 0.5 convention on both sides)."""
    p = os.path.join(tempfile.mkdtemp(prefix="qp_tiny_"), "t.mps")
    open(p, "w").write("NAME t\nROWS\n N obj\nCOLUMNS\n    x obj 1\nBOUNDS\n LO bnd x 1\nQUADOBJ\n    x x 2\nENDATA\n")
    h = run_highs(p)
    js = p + ".json"
    subprocess.run([engine, p, "--json", js], capture_output=True, timeout=60)
    t = json.load(open(js))["objective"]
    assert close(h[1], 2) and close(t, 2), "0.5 convention mismatch: HiGHS %r, taral %r" % (h, t)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--engine", default=os.path.join(os.path.dirname(__file__), "..", "out", "taral"))
    ap.add_argument("--per-suite", type=int, default=25)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--keep", help="directory to keep MPS files of non-ok cases")
    a = ap.parse_args()
    engine = os.path.abspath(a.engine)
    tiny(engine)
    tmp = tempfile.mkdtemp(prefix="qp_tests_")
    jobs = []
    for s, (name, gen) in enumerate(SUITES.items()):
        rng = np.random.default_rng([a.seed, s])
        for k in range(a.per_suite):
            jobs.append((name, gen(rng, k)))

    def one(job):
        suite, case = job
        path = os.path.join(tmp, "%s_%s.mps" % (suite, case.name))
        write_mps(case, path)
        hst, hobj = run_highs(path)
        r, x = run_taral(engine, path)
        verdict, why = judge(case, r, x, hst, hobj)
        return suite, case, verdict, why, r, hst, hobj, path

    tally, wrong = {}, 0
    with ThreadPoolExecutor(a.workers) as ex:
        for suite, case, verdict, why, r, hst, hobj, path in ex.map(one, jobs):
            t = tally.setdefault(suite, {"ok": 0, "undecided": 0, "wrong": 0, "status": {}})
            t[verdict] += 1
            t["status"][r["status"]] = t["status"].get(r["status"], 0) + 1
            if verdict != "ok":
                print("%-9s %-12s %-14s taral=%s obj=%s | HiGHS %s %s | %s" % (
                    verdict.upper(), suite, case.name, r["status"], r.get("objective"), hst, hobj, why))
                if a.keep:
                    os.makedirs(a.keep, exist_ok=True)
                    os.replace(path, os.path.join(a.keep, os.path.basename(path)))
            wrong += verdict == "wrong"
    total = {"ok": 0, "undecided": 0, "wrong": 0}
    for suite, t in tally.items():
        print("%-13s ok %3d  undecided %3d  wrong %3d   statuses %s" % (
            suite, t["ok"], t["undecided"], t["wrong"], dict(sorted(t["status"].items()))))
        for k in total:
            total[k] += t[k]
    print("TOTAL %d cases: ok %d, undecided %d, wrong %d (seed %d)" % (
        sum(total.values()), total["ok"], total["undecided"], total["wrong"], a.seed))
    sys.exit(1 if wrong else 0)


if __name__ == "__main__":
    main()
