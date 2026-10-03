#!/usr/bin/env python3
"""Large-magnitude MILP generator (supplementary check, not an acceptance set).

Takes the integer-feasible random models of benchmarks/milp_stress.py's shape (5-25 columns, 3-15 rows, a known
integer point) and rescales them so that coefficients, right-hand sides, costs or variable values are large:
  obj    cost vector x 10^3..10^9           row    each row x 10^0..10^7          cost   per-column cost x 10^0..10^8
  bound  columns shifted by up to 10^6      mid    rows x 10^0..10^4, cost x 10^0..10^5
  wide   rows x 10^-2..10^5, cost x 10^-1..10^4 per column
  mixed  rows x 10^-3..10^7, cost x 10^-2..10^8 (root LP failures live here, see NOTES.md)
  all    obj + row together                 range  integer columns with ranges up to 10^6 around +-10^5
Each model is solved by the engine and by HiGHS (mip gaps 1e-9) on the same MPS. Only the interesting cases are
printed: status differences, objective differences, and any run that ends numerical_failure or with unresolved
nodes. Seed 504, case 17 of the `all` family is the model in cpp-engine/tests/milp/big_values_unresolved.mps.

  python benchmarks/milp_bigvals.py --engine OUT/taral --family all --seeds 504 505 --n 40 [--keep DIR]
"""
import argparse, json, os, random, subprocess, sys, tempfile
import numpy as np
import highspy


def gen(rng, big_range):
    n = rng.randint(5, 25); m = rng.randint(3, 15)
    if big_range:
        lo = [rng.randint(-10**5, 10**5) for _ in range(n)]; up = [l + rng.randint(1, 10 ** rng.randint(1, 6)) for l in lo]
    else:
        lo = [rng.randint(-5, 2) for _ in range(n)]; up = [l + rng.randint(1, 15) for l in lo]
    ints = [rng.random() < 0.7 for _ in range(n)]
    pt = [rng.randint(l, u) for l, u in zip(lo, up)]
    A = np.zeros((m, n)); typ = []; rhs = []
    for i in range(m):
        for j in range(n):
            if rng.random() < 0.5: A[i, j] = rng.randint(-5, 9)
        act = A[i] @ pt; t = rng.choice("LLGE"); typ.append(t)
        rhs.append(act + rng.randint(0, 6) if t == "L" else act - rng.randint(0, 6) if t == "G" else act)
    c = np.array([rng.randint(-9, 9) for _ in range(n)], float)
    return n, m, np.array(lo, float), np.array(up, float), ints, A, typ, np.array(rhs, float), c, rng.random() < 0.5


def scale(rng, kind, n, m, lo, up, A, rhs, c):
    if kind in ("obj", "all"): c = c * 10.0 ** rng.randint(3, 9)
    if kind in ("row", "all"):
        s = 10.0 ** np.array([rng.randint(0, 7) for _ in range(m)]); A = A * s[:, None]; rhs = rhs * s
    if kind == "mid":
        s = 10.0 ** np.array([rng.randint(0, 4) for _ in range(m)]); A = A * s[:, None]; rhs = rhs * s
        c = c * 10.0 ** rng.randint(0, 5)
    if kind == "wide":
        s = 10.0 ** np.array([rng.randint(-2, 5) for _ in range(m)]); A = A * s[:, None]; rhs = rhs * s
        c = c * 10.0 ** np.array([rng.randint(-1, 4) for _ in range(n)])
    if kind == "mixed":
        s = 10.0 ** np.array([rng.randint(-3, 7) for _ in range(m)]); A = A * s[:, None]; rhs = rhs * s
        c = c * 10.0 ** rng.randint(-2, 8)
    if kind == "cost": c = c * 10.0 ** np.array([rng.randint(0, 8) for _ in range(n)])
    if kind == "bound":
        sh = np.array([rng.randint(0, 10**6) for _ in range(n)], float); rhs = rhs + A @ sh; lo = lo + sh; up = up + sh
    return A, rhs, lo, up, c


def write(path, n, m, lo, up, ints, A, typ, rhs, c, maxi):
    L = ["NAME B"] + (["OBJSENSE", "    MAX"] if maxi else []) + ["ROWS", " N OBJ"] + [f" {t} R{i}" for i, t in enumerate(typ)] + ["COLUMNS"]
    inint = False
    for j in range(n):
        if ints[j] and not inint: L.append("    MARKER 'MARKER' 'INTORG'"); inint = True
        if not ints[j] and inint: L.append("    MARKER 'MARKER' 'INTEND'"); inint = False
        if c[j] != 0: L.append(f"    X{j} OBJ {float(c[j])!r}")
        for i in range(m):
            if A[i, j] != 0: L.append(f"    X{j} R{i} {float(A[i, j])!r}")
    if inint: L.append("    MARKER 'MARKER' 'INTEND'")
    L.append("RHS"); L += [f"    RHS R{i} {float(rhs[i])!r}" for i in range(m)]
    L.append("BOUNDS")
    for j in range(n): L += [f" LO BND X{j} {float(lo[j])!r}", f" UP BND X{j} {float(up[j])!r}"]
    L.append("ENDATA"); open(path, "w").write("\n".join(L) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--engine", required=True)
    ap.add_argument("--family", required=True, choices="obj row cost bound mid wide mixed all range".split())
    ap.add_argument("--seeds", type=int, nargs="+", required=True)
    ap.add_argument("--n", type=int, default=100)
    ap.add_argument("--time-limit", type=float, default=10)
    ap.add_argument("--keep", default=None, help="directory for the MPS files (default: a temporary one)")
    a = ap.parse_args()
    d = a.keep or tempfile.mkdtemp(); os.makedirs(d, exist_ok=True)
    stats = {}; flagged = 0
    for seed in a.seeds:
        rng = random.Random(seed)
        for k in range(a.n):
            n, m, lo, up, ints, A, typ, rhs, c, maxi = gen(rng, a.family == "range")
            A, rhs, lo, up, c = scale(rng, a.family, n, m, lo, up, A, rhs, c)
            p = os.path.join(d, f"s{seed}_{k}.mps"); write(p, n, m, lo, up, ints, A, typ, rhs, c, maxi)
            h = highspy.Highs(); h.setOptionValue("output_flag", False); h.setOptionValue("time_limit", 30.0)
            h.setOptionValue("mip_rel_gap", 1e-9); h.setOptionValue("mip_abs_gap", 1e-9); h.readModel(p); h.run()
            hs = h.modelStatusToString(h.getModelStatus()); ho = h.getInfo().objective_function_value
            js = p + ".json"
            subprocess.run([a.engine, p, "--time-limit", str(a.time_limit), "--json", js, "--sol", p + ".sol"], capture_output=True)
            r = json.load(open(js)); st = r["status"]; stats[st] = stats.get(st, 0) + 1
            ref = "optimal" if hs == "Optimal" else "infeasible" if hs == "Infeasible" else hs
            tag = "ok"
            if st != ref: tag = "STATUS"
            elif st == "optimal" and abs(r["objective"] - ho) > 1e-6 * max(1, abs(ho)): tag = "OBJECTIVE"
            if tag != "ok" or st == "numerical_failure" or (r.get("unresolved_nodes") or 0) > 0:
                flagged += 1
                print(f"{tag} seed={seed} case={k} engine={st} {r.get('objective')} nodes={r.get('nodes')} "
                      f"unresolved={r.get('unresolved_nodes')} reference={hs} {ho} msg={r.get('message')!r} file={p}", flush=True)
    print(f"DONE family={a.family} statuses={stats} flagged={flagged}")


main()
