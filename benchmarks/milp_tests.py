#!/usr/bin/env python3
"""Seeded randomized MILP suites for the C++ engine, checked against HiGHS on the same MPS file.

Every case is written to MPS, solved by `taral` (CLI, JSON + .sol) and by HiGHS (highspy, tight gaps).
Independently of HiGHS, every reported solution is re-checked here against the generated data
(rows, bounds, integrality, objective). Suites:
  a mixed        mixed continuous + integer, L/G/E/ranged rows
  b binary       knapsack (single and multi constraint) and assignment with a side constraint
  c infeasible   integer-infeasible models whose LP relaxation is feasible (parity / narrow bounds)
  d unbounded    unbounded MILPs, and unbounded relaxations with no integer point
  e weak         set covering and big-M fixed-charge facility location
  f maximize     OBJSENSE MAX versions of mixed models
  g general      general integers with negative / fractional / one-sided bounds
  h limits       node and time limits: incumbent >= optimum >= reported bound (minimisation sense)

A case is WRONG when taral claims something false (an optimum that differs from HiGHS, a status that
contradicts the known truth, an infeasible point, an invalid bound). A case is UNDECIDED when taral
honestly reports it could not decide (limit / unbounded_relaxation); those are counted separately.

Usage: ~/.venvs/taral-gpu/bin/python benchmarks/milp_tests.py [--engine out/taral] [--per-suite 40] [--seed 1]
Exit code 1 on any WRONG case.
"""
import argparse, json, math, os, subprocess, sys, tempfile
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import highspy

INF = math.inf


class Case:
    def __init__(self, name, c, A, rlo, rup, lo, up, isint, maximize=False, truth=None, args=()):
        self.name, self.c, self.A = name, np.asarray(c, float), np.asarray(A, float)
        self.rlo, self.rup = np.asarray(rlo, float), np.asarray(rup, float)
        self.lo, self.up = np.asarray(lo, float), np.asarray(up, float)
        self.isint, self.maximize = np.asarray(isint, bool), maximize
        self.truth = truth  # known status by construction, or None (HiGHS decides)
        self.args = list(args)  # extra CLI arguments (limits)


def num(v):
    return repr(float(v)) if v != int(v) else str(int(v))


def write_mps(case, path):
    m, n = case.A.shape
    out = ["NAME " + case.name]
    if case.maximize:
        out += ["OBJSENSE", "    MAX"]
    out.append("ROWS")
    out.append(" N obj")
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
    in_int = False
    for j in range(n):
        if case.isint[j] != in_int:
            out.append("    M%d 'MARKER' '%s'" % (j, "INTORG" if case.isint[j] else "INTEND"))
            in_int = case.isint[j]
        out.append("    c%d obj %s" % (j, num(case.c[j])))
        for i in range(m):
            if case.A[i, j] != 0:
                out.append("    c%d r%d %s" % (j, i, num(case.A[i, j])))
    if in_int:
        out.append("    MEND 'MARKER' 'INTEND'")
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
    if case.isint.any():
        w = max(w, float(np.max(np.abs(x[case.isint] - np.round(x[case.isint])))))
    return w


def run_highs(path):
    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    h.setOptionValue("mip_rel_gap", 1e-9)
    h.setOptionValue("mip_abs_gap", 1e-9)
    h.setOptionValue("time_limit", 120.0)
    h.setOptionValue("threads", 1)
    h.readModel(path)
    h.run()
    st = h.modelStatusToString(h.getModelStatus())
    info = h.getInfo()
    return st, info.objective_function_value


def run_taral(engine, path, args):
    js, sol = path + ".json", path + ".sol"
    for p in (js, sol):
        if os.path.exists(p):
            os.remove(p)
    subprocess.run([engine, path, "--time-limit", "60", "--json", js, "--sol", sol] + args,
                   capture_output=True, timeout=600)
    with open(js) as f:
        r = json.load(f)
    x = None
    if os.path.exists(sol):
        vals = dict(line.split() for line in open(sol))
        x = np.array([float(vals["c%d" % j]) for j in range(len(vals))])
    return r, x


def close(a, b):
    return abs(a - b) <= 1e-6 * max(1.0, abs(a), abs(b))


def judge(case, r, x, hst, hobj):
    """Return (verdict, reason): verdict in ok / undecided / wrong."""
    st = r["status"]
    sense = -1 if case.maximize else 1
    truth = case.truth or {"Optimal": "optimal", "Infeasible": "infeasible", "Unbounded": "unbounded"}.get(hst)
    if r["has_solution"]:
        if x is None:
            return "wrong", "has_solution without a .sol file"
        v = violation(case, x)
        if v > 1e-6:
            return "wrong", "reported point violates the model by %.3g" % v
        if st != "unbounded" and r["objective"] is not None and not close(r["objective"], float(case.c @ x)):
            return "wrong", "objective %r does not match the point (%r)" % (r["objective"], float(case.c @ x))
    if st == "optimal":
        if truth != "optimal":
            return "wrong", "optimal but truth/HiGHS says %s" % (truth or hst)
        if not r["has_solution"] or not close(r["objective"], hobj):
            return "wrong", "objective %r vs HiGHS %r" % (r["objective"], hobj)
        if r["gap"] is None or r["gap"] > 1e-6 + 1e-12:
            return "wrong", "optimal with gap %r" % r["gap"]
        return "ok", ""
    if st in ("infeasible", "unbounded"):
        if truth == st:
            return "ok", ""
        if st == "unbounded" and hst == "Primal infeasible or unbounded" and r["has_solution"]:
            return "ok", "HiGHS says infeasible-or-unbounded; feasible point shown"
        return "wrong", "%s but truth/HiGHS says %s" % (st, truth or hst)
    if st in ("time_limit", "node_limit", "unbounded_relaxation", "numerical_failure"):
        if truth == "infeasible" and r["has_solution"]:
            return "wrong", "point reported for an infeasible model"
        if truth == "optimal" and st != "unbounded_relaxation":
            b = r["best_bound"]
            if b is not None and sense * b > sense * hobj + 1e-6 * max(1, abs(hobj)):
                return "wrong", "bound %r cuts off the optimum %r" % (b, hobj)
            if r["has_solution"] and sense * r["objective"] < sense * hobj - 1e-6 * max(1, abs(hobj)):
                return "wrong", "incumbent %r better than the optimum %r" % (r["objective"], hobj)
        if truth == "optimal" and st == "unbounded_relaxation":
            return "wrong", "unbounded_relaxation on a model with a finite optimum"
        return "undecided", st
    return "wrong", "unknown status " + st


# ---------------------------------------------------------------- generators
def g_mixed(rng, k, maximize=False):
    n, m = rng.integers(6, 18), rng.integers(4, 12)
    A = rng.integers(-6, 7, (m, n)) * (rng.random((m, n)) < 0.6)
    isint = rng.random(n) < 0.6
    lo = np.zeros(n); up = rng.integers(3, 15, n).astype(float)
    x0 = np.where(isint, rng.integers(0, 3, n), rng.random(n) * 3)  # a feasible point
    act = A @ x0
    rlo, rup = np.full(m, -INF), np.full(m, INF)
    for i in range(m):
        t = rng.integers(0, 10)
        if t < 4: rup[i] = math.floor(act[i]) + rng.integers(0, 6)
        elif t < 7: rlo[i] = math.ceil(act[i]) - rng.integers(0, 6)
        elif t < 8: rlo[i] = rup[i] = act[i]  # equality through the feasible point
        else: rlo[i], rup[i] = math.ceil(act[i]) - rng.integers(0, 4), math.floor(act[i]) + rng.integers(1, 5)
    rlo = np.minimum(rlo, act); rup = np.maximum(rup, act)
    c = rng.integers(-10, 11, n).astype(float)
    return Case("mixed%d" % k, c, A, rlo, rup, lo, up, isint, maximize)


def g_binary(rng, k):
    if k % 2 == 0:  # knapsack, 1-3 constraints, maximise value
        n, m = rng.integers(10, 30), rng.integers(1, 4)
        w = rng.integers(5, 60, (m, n))
        cap = np.floor(w.sum(1) * rng.uniform(0.25, 0.6))
        v = rng.integers(5, 80, n).astype(float)
        return Case("knap%d" % k, v, w, np.full(m, -INF), cap, np.zeros(n), np.ones(n), np.ones(n, bool), True)
    p = rng.integers(3, 7)  # assignment p x p with a budget side constraint
    n = p * p
    cost = rng.integers(1, 30, n).astype(float)
    rows, rlo, rup = [], [], []
    for i in range(p):
        r = np.zeros(n); r[i * p:(i + 1) * p] = 1; rows.append(r); rlo.append(1); rup.append(1)
        r = np.zeros(n); r[i::p] = 1; rows.append(r); rlo.append(1); rup.append(1)
    t = rng.integers(1, 20, n)
    rows.append(t); rlo.append(-INF); rup.append(float(np.sort(t)[:p].sum() + rng.integers(2 * p, 8 * p)))
    return Case("assign%d" % k, cost, np.array(rows), rlo, rup, np.zeros(n), np.ones(n), np.ones(n, bool))


def g_infeasible(rng, k):
    base = g_mixed(rng, k)
    n = len(base.c)
    if k % 2 == 0:  # parity: sum of even multiples of bounded integers equals an odd number
        j = [n, n + 1, n + 2]
        A = np.hstack([base.A, np.zeros((base.A.shape[0], 3))])
        row = np.zeros(n + 3); row[j] = 2 * rng.integers(1, 4, 3)
        A = np.vstack([A, row])
        odd = 2 * rng.integers(2, 6) + 1
        return Case("parity%d" % k, np.append(base.c, [1, 1, 1]), A, np.append(base.rlo, odd), np.append(base.rup, odd),
                    np.append(base.lo, [0, 0, 0]), np.append(base.up, [9, 9, 9]), np.append(base.isint, [True] * 3),
                    truth="infeasible")
    # narrow window: a*x + a*y in [a*q + 1, a*q + a - 1] has no integer solution, LP feasible
    a, q = int(rng.integers(3, 8)), int(rng.integers(1, 5))
    A = np.hstack([base.A, np.zeros((base.A.shape[0], 2))])
    row = np.zeros(n + 2); row[[n, n + 1]] = a
    A = np.vstack([A, row])
    return Case("window%d" % k, np.append(base.c, [1, -1]), A, np.append(base.rlo, a * q + 1),
                np.append(base.rup, a * q + a - 1), np.append(base.lo, [0, 0]), np.append(base.up, [10, 10]),
                np.append(base.isint, [True, True]), truth="infeasible")


def g_unbounded(rng, k):
    base = g_mixed(rng, k)
    n, m = len(base.c), base.A.shape[0]
    if k % 2 == 0:  # an integer column that only loosens <= rows, improving cost, no upper bound
        col = np.where(np.isinf(base.rlo), -rng.integers(0, 3, m), 0)
        A = np.hstack([base.A, col[:, None]])
        return Case("unb%d" % k, np.append(base.c, -1), A, base.rlo, base.rup, np.append(base.lo, 0),
                    np.append(base.up, INF), np.append(base.isint, True), truth="unbounded")
    # relaxation unbounded through a free continuous column, integers infeasible by parity (bounded)
    A = np.hstack([base.A, np.zeros((m, 3))])
    row = np.zeros(n + 3); row[[n, n + 1]] = 2
    A = np.vstack([A, row])
    return Case("unbrel%d" % k, np.append(base.c, [0, 0, -1]), A, np.append(base.rlo, 7), np.append(base.rup, 7),
                np.append(base.lo, [0, 0, -INF]), np.append(base.up, [6, 6, INF]),
                np.append(base.isint, [True, True, False]), truth="infeasible")


def g_weak(rng, k):
    if k % 2 == 0:  # set covering
        n, m = rng.integers(15, 40), rng.integers(10, 30)
        A = (rng.random((m, n)) < 0.15).astype(float)
        for i in range(m):
            if not A[i].any(): A[i, rng.integers(0, n)] = 1
        c = rng.integers(1, 20, n).astype(float)
        return Case("cover%d" % k, c, A, np.ones(m), np.full(m, INF), np.zeros(n), np.ones(n), np.ones(n, bool))
    f, d = rng.integers(3, 7), rng.integers(4, 10)  # facility location, aggregated big-M
    dem = rng.integers(5, 30, d).astype(float)
    M = dem.sum()
    n = f + f * d
    c = np.concatenate([rng.integers(50, 300, f), rng.integers(1, 20, f * d)]).astype(float)
    rows, rlo, rup = [], [], []
    for j in range(d):  # demand
        r = np.zeros(n); r[f + np.arange(f) * d + j] = 1; rows.append(r); rlo.append(dem[j]); rup.append(INF)
    for i in range(f):  # sum_j x_ij <= M y_i
        r = np.zeros(n); r[f + i * d:f + (i + 1) * d] = 1; r[i] = -M; rows.append(r); rlo.append(-INF); rup.append(0)
    lo, up = np.zeros(n), np.concatenate([np.ones(f), np.full(f * d, INF)])
    return Case("facil%d" % k, c, np.array(rows), rlo, rup, lo, up, np.arange(n) < f)


def g_maximize(rng, k):
    return g_mixed(rng, k, maximize=True)


def g_general(rng, k):
    n, m = rng.integers(5, 14), rng.integers(3, 9)
    A = rng.integers(-5, 6, (m, n)) * (rng.random((m, n)) < 0.7)
    isint = rng.random(n) < 0.8
    lo = rng.integers(-12, 0, n).astype(float) + np.where(rng.random(n) < 0.3, -0.5, 0)
    up = lo + rng.integers(2, 20, n) + np.where(rng.random(n) < 0.3, 0.3, 0)
    lo[rng.random(n) < 0.15] = -INF  # one-sided; rows below keep the box of x0 feasible
    x0 = np.where(isint, np.ceil(np.maximum(lo, -5)), np.maximum(lo, -5) + 0.3)
    x0 = np.minimum(x0, np.floor(up))
    act = A @ x0
    rlo = act - rng.integers(0, 8, m); rup = act + rng.integers(0, 8, m)
    rlo[rng.random(m) < 0.3] = -INF
    # keep the free-below columns bounded through a box row so the optimum is finite
    rows = [A]; blo = [rlo]; bup = [rup]
    for j in np.where(lo == -INF)[0]:
        r = np.zeros(n); r[j] = 1; rows.append(r[None]); blo.append([x0[j] - 7]); bup.append([INF])
    c = rng.integers(-9, 10, n).astype(float)
    return Case("gen%d" % k, c, np.vstack(rows), np.concatenate(blo), np.concatenate(bup), lo, up, isint,
                maximize=bool(k % 3 == 0))


def g_limits(rng, k):
    if k % 2 == 0:
        base = g_weak(rng, 0) if k % 4 else g_binary(rng, 0)
        base.args = ["--node-limit", str(int(rng.integers(1, 8)))]
    else:  # larger covering model under a short time limit
        n, m = 120, 80
        A = (rng.random((m, n)) < 0.06).astype(float)
        for i in range(m):
            A[i, rng.integers(0, n)] = 1
        base = Case("", rng.integers(1, 40, n).astype(float), A, np.ones(m), np.full(m, INF), np.zeros(n), np.ones(n),
                    np.ones(n, bool))
        base.args = ["--time-limit", "0.01"]
    base.name = "lim%d" % k
    return base


SUITES = {"a_mixed": g_mixed, "b_binary": g_binary, "c_infeasible": g_infeasible, "d_unbounded": g_unbounded,
          "e_weak": g_weak, "f_maximize": g_maximize, "g_general": g_general, "h_limits": g_limits}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--engine", default=os.path.join(os.path.dirname(__file__), "..", "out", "taral"))
    ap.add_argument("--per-suite", type=int, default=40)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--keep", help="directory to keep MPS files of wrong cases")
    a = ap.parse_args()
    engine = os.path.abspath(a.engine)
    tmp = tempfile.mkdtemp(prefix="milp_tests_")
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
        r, x = run_taral(engine, path, case.args)
        verdict, why = judge(case, r, x, hst, hobj)
        return suite, case, verdict, why, r, hst, hobj, path

    tally = {}
    wrong = 0
    with ThreadPoolExecutor(a.workers) as ex:
        for suite, case, verdict, why, r, hst, hobj, path in ex.map(one, jobs):
            t = tally.setdefault(suite, {"ok": 0, "undecided": 0, "wrong": 0, "status": {}})
            t[verdict] += 1
            t["status"][r["status"]] = t["status"].get(r["status"], 0) + 1
            if verdict != "ok":
                print("%-9s %-12s %-14s taral=%s obj=%s bound=%s nodes=%s | HiGHS %s %s | %s" % (
                    verdict.upper(), suite, case.name, r["status"], r["objective"], r["best_bound"], r["nodes"],
                    hst, hobj, why))
            if verdict == "wrong":
                wrong += 1
                if a.keep:
                    os.makedirs(a.keep, exist_ok=True)
                    os.replace(path, os.path.join(a.keep, os.path.basename(path)))
    total = {"ok": 0, "undecided": 0, "wrong": 0}
    for suite, t in tally.items():
        print("%-13s ok %3d  undecided %3d  wrong %3d   statuses %s" % (suite, t["ok"], t["undecided"], t["wrong"],
                                                                     dict(sorted(t["status"].items()))))
        for k in total:
            total[k] += t[k]
    print("TOTAL %d cases: ok %d, undecided %d, wrong %d (local run, seed %d)" % (
        sum(total.values()), total["ok"], total["undecided"], total["wrong"], a.seed))
    sys.exit(1 if wrong else 0)


if __name__ == "__main__":
    main()
