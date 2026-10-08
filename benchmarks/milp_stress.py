#!/usr/bin/env python3
"""Supplementary MILP stress check (not an acceptance set): larger random models than milp_check.py, with
equality rows built around a known integer point, negative lower bounds, general integers, and big-M style
rows. Graded against HiGHS (scipy.optimize.milp on the generator data) and an independent check of the
engine's point. Prints one line; exit 1 on any wrong answer.

  python benchmarks/milp_stress.py --engine OUT/taral --seeds 1 2 3 --n 100   (add --audit-prop to the engine command to cross-check propagation prunes)
"""
import argparse, json, os, random, subprocess, sys, tempfile
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tools"))
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tools", "milp_check.py")).read()
ns = {"__name__": "milp_check_helpers"}
exec(compile(src.rsplit("\nmain()", 1)[0], "milp_check.py", "exec"), ns)
write_mps = ns["write_mps"]


def gen(rng, k):
    n = rng.randint(5, 25); m = rng.randint(3, 15); cols = [f"X{j}" for j in range(n)]
    lo = {c: rng.randint(-5, 2) for c in cols}; bounds = {c: (lo[c], lo[c] + rng.randint(1, 15)) for c in cols}
    ints = [rng.random() < 0.7 for _ in cols]
    pt = {c: rng.randint(bounds[c][0], bounds[c][1]) for c in cols}  # a feasible integer point
    rows, rhs = [], {}
    for i in range(m):
        row = {c: rng.randint(-5, 9) for c in cols if rng.random() < 0.5}
        act = sum(v * pt[c] for c, v in row.items()); r = f"R{i}"; t = rng.choice("LLGE")
        rows.append((r, t, row))
        rhs[r] = act + rng.randint(0, 6) if t == "L" else act - rng.randint(0, 6) if t == "G" else act
        if rng.random() < 0.15:  # sometimes cut the known point off
            rhs[r] += rng.choice([-1, 1]) * rng.randint(3, 30)
    obj = [rng.randint(-9, 9) for _ in cols]
    return dict(name=f"STR{k}", rows=rows, cols=cols, obj=obj, rhs=rhs, bounds=bounds, ints=ints, maxi=rng.random() < 0.5)


def reference(g):
    cols = g["cols"]; n = len(cols); c = np.array(g["obj"], float) * (-1 if g["maxi"] else 1)
    A = np.zeros((len(g["rows"]), n)); lo = np.zeros(len(g["rows"])); up = np.zeros(len(g["rows"]))
    for i, (r, t, row) in enumerate(g["rows"]):
        for j, cn in enumerate(cols): A[i, j] = row.get(cn, 0)
        b = g["rhs"][r]; lo[i], up[i] = (-np.inf, b) if t == "L" else (b, np.inf) if t == "G" else (b, b)
    res = milp(c, constraints=LinearConstraint(A, lo, up),
               bounds=Bounds([g["bounds"][x][0] for x in cols], [g["bounds"][x][1] for x in cols]),
               integrality=np.array([1 if x else 0 for x in g["ints"]]))
    if res.status == 0: return "optimal", float(res.fun) * (-1 if g["maxi"] else 1)
    if res.status == 2: return "infeasible", 0.0
    return "other", 0.0


def verify(g, solpath, tol=1e-6):
    x = {}
    for line in open(solpath):
        nm, v = line.split(); x[nm] = float(v)
    for j, cn in enumerate(g["cols"]):
        lo, up = g["bounds"][cn]; v = x.get(cn, 0.0)
        if v < lo - tol or v > up + tol or (g["ints"][j] and abs(v - round(v)) > tol): return None
    for r, t, row in g["rows"]:
        a = sum(v * x.get(cn, 0.0) for cn, v in row.items()); b = g["rhs"][r]
        if (t in "LE" and a > b + tol) or (t in "GE" and a < b - tol): return None
    return sum(o * x.get(cn, 0.0) for o, cn in zip(g["obj"], g["cols"]))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--engine", required=True); ap.add_argument("--seeds", type=int, nargs="+", required=True)
    ap.add_argument("--n", type=int, default=100); ap.add_argument("--time-limit", type=float, default=30)
    a = ap.parse_args(); d = tempfile.mkdtemp(); bad = nodes = cases = unsolved = pruned = 0; st = {}; audit = [0] * 10
    for seed in a.seeds:
        rng = random.Random(seed)
        for k in range(a.n):
            g = gen(rng, k); p = os.path.join(d, "s%d_%s.mps" % (seed, g["name"])); write_mps(p, **g)
            hs, ho = reference(g); js, sp = p + ".json", p + ".sol"
            subprocess.run([a.engine, p, "--time-limit", str(a.time_limit), "--json", js, "--sol", sp], capture_output=True)
            r = json.load(open(js)); cases += 1; nodes += r.get("nodes") or 0; pruned += r.get("prop_crossed") or 0
            audit = [x + y for x, y in zip(audit, r.get("audit") or [0] * 10)]
            st[r["status"]] = st.get(r["status"], 0) + 1
            if r["status"] in ("time_limit", "node_limit"): unsolved += 1; continue
            ok = r["status"] == hs
            if ok and hs == "optimal":
                o = verify(g, sp); ok = o is not None and abs(o - ho) <= 1e-6 * max(1, abs(ho)) and abs(o - r["objective"]) <= 1e-6 * max(1, abs(o))
            if not ok:
                bad += 1; print("MISMATCH seed", seed, g["name"], "reference", hs, ho, "engine", r["status"], r.get("objective"))
    print(f"STRESS cases={cases} wrong={bad} unsolved={unsolved} total_nodes={nodes} prop_crossed={pruned} audit[nodes,lp_inf,lp_feas,lp_other,int_empty,int_feasible,undecided,cutoff_only,rc_checked,rc_bad]={audit} statuses={st}"); sys.exit(1 if bad else 0)


main()
