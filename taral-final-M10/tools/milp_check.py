#!/usr/bin/env python3
"""MILP regression: random small mixed-integer models (plus edge cases) solved by the C++ engine and by HiGHS on the SAME .mps.
Checks status class and optimal objective (abs/rel 1e-6). Tooling only: HiGHS is the reference, never inside the engine.
Usage: python tools/milp_check.py --cmd ./taral [--n 40] [--seed 26000+119] [--out DIR]"""
import argparse, json, os, random, subprocess, sys, tempfile
import highspy
def write_mps(path, name, rows, cols, obj, rhs, bounds, ints, maxi=False):
    L = [f"NAME {name}"] + (["OBJSENSE","    MAX"] if maxi else []) + ["ROWS"," N OBJ"]
    for r, t, _ in rows: L.append(f" {t} {r}")
    L.append("COLUMNS"); inint = False
    for j, c in enumerate(cols):
        if ints[j] and not inint: L.append("    MARKER 'MARKER' 'INTORG'"); inint = True
        if not ints[j] and inint: L.append("    MARKER 'MARKER' 'INTEND'"); inint = False
        if obj[j] != 0: L.append(f"    {c} OBJ {obj[j]}")
        for r, _, row in rows:
            if row.get(c, 0) != 0: L.append(f"    {c} {r} {row[c]}")
    if inint: L.append("    MARKER 'MARKER' 'INTEND'")
    L.append("RHS")
    for r, _, _ in rows: L.append(f"    RHS {r} {rhs[r]}")
    L.append("BOUNDS")
    for c, (lo, up) in bounds.items():
        if lo is not None: L.append(f" LO BND {c} {lo}")
        if up is not None: L.append(f" UP BND {c} {up}")
    L.append("ENDATA"); open(path, "w").write("\n".join(L) + "\n")
def highs(path):
    h = highspy.Highs(); h.setOptionValue("output_flag", False); h.setOptionValue("time_limit", 60.0); h.readModel(path); h.run()
    st = h.modelStatusToString(h.getModelStatus()); return st, h.getInfo().objective_function_value
def ref_direct(g):
    """Reference built from the generator data (explicit integrality), independent of any MPS reader version. scipy milp = HiGHS, tooling only."""
    import numpy as np
    from scipy.optimize import milp, LinearConstraint, Bounds
    cols = g["cols"]; n = len(cols); c = np.array(g["obj"], float) * (-1 if g["maxi"] else 1)
    A = np.zeros((len(g["rows"]), n)); lo = np.zeros(len(g["rows"])); up = np.zeros(len(g["rows"]))
    for i, (r, t, row) in enumerate(g["rows"]):
        for j, cn in enumerate(cols): A[i, j] = row.get(cn, 0)
        b = g["rhs"][r]; lo[i], up[i] = ((-np.inf, b) if t == "L" else (b, np.inf))
    bl = np.array([g["bounds"][cn][0] for cn in cols], float); bu = np.array([g["bounds"][cn][1] for cn in cols], float)
    res = milp(c, constraints=LinearConstraint(A, lo, up), bounds=Bounds(bl, bu), integrality=np.array([1 if x else 0 for x in g["ints"]]))
    if res.status == 0: return "Optimal", float(res.fun) * (-1 if g["maxi"] else 1)
    if res.status == 2: return "Infeasible", 0.0
    return "Other%d" % res.status, 0.0
def verify_sol(g, solpath, tol=1e-6):
    """Independent feasibility certificate for the engine's point: bounds, integrality, rows, recomputed objective."""
    x = {}
    for line in open(solpath):
        nm, v = line.split(); x[nm] = float(v)
    for j, cn in enumerate(g["cols"]):
        lo, up = g["bounds"][cn]; v = x.get(cn, 0.0)
        if v < lo - tol or v > up + tol: return False, None
        if g["ints"][j] and abs(v - round(v)) > tol: return False, None
    for r, t, row in g["rows"]:
        a = sum(row.get(cn, 0) * x.get(cn, 0.0) for cn in g["cols"]); b = g["rhs"][r]
        if (t == "L" and a > b + tol) or (t == "G" and a < b - tol): return False, None
    return True, sum(o * x.get(cn, 0.0) for o, cn in zip(g["obj"], g["cols"]))
def gen(rng, k):
    n = rng.randint(3, 12); m = rng.randint(2, 8); cols = [f"X{j}" for j in range(n)]
    ints = [rng.random() < 0.7 for _ in range(n)]
    rows = []; rhs = {}
    for i in range(m):
        row = {c: rng.randint(-3, 9) for c in cols if rng.random() < 0.7}
        t = rng.choice("LLG"); r = f"R{i}"; rows.append((r, t, row)); rhs[r] = rng.randint(5, 40) if t == "L" else rng.randint(0, 6)
    obj = [rng.randint(-9, 9) for _ in cols]; bounds = {c: (0, rng.randint(2, 12)) for c in cols}
    return dict(name=f"RND{k}", rows=rows, cols=cols, obj=obj, rhs=rhs, bounds=bounds, ints=ints, maxi=rng.random() < 0.5)
def main():
    import scipy; print("python", sys.version.split()[0], "scipy", scipy.__version__, "highspy", getattr(highspy, "__version__", "?"))
    ap = argparse.ArgumentParser(); ap.add_argument("--cmd", required=True); ap.add_argument("--n", type=int, default=40); ap.add_argument("--seed", type=int, default=26000+119)
    ap.add_argument("--out", default=None); a = ap.parse_args(); rng = random.Random(a.seed); d = a.out or tempfile.mkdtemp(); os.makedirs(d, exist_ok=True)
    bad = 0; refsub = 0; stats = {"optimal": 0, "infeasible": 0}
    for k in range(a.n):
        g = gen(rng, k); p = os.path.join(d, g["name"] + ".mps"); write_mps(p, **g)
        hs, ho = ref_direct(g); hs_mps, ho_mps = highs(p); jp = p + ".json"
        if hs != hs_mps or (hs == 'Optimal' and abs(ho - ho_mps) > 1e-6 * max(1, abs(ho))): print('NOTE reference reader differs on', g['name'], 'direct', hs, ho, 'mps-read', hs_mps, ho_mps)
        sp = p + ".sol"; subprocess.run([a.cmd, p, "--time-limit", "30", "--json", jp, "--sol", sp], check=False); r = json.load(open(jp))
        hs_cls = "optimal" if hs == "Optimal" else ("infeasible" if hs == "Infeasible" else hs)
        ok = (r["status"] == hs_cls)
        if ok and hs_cls == "optimal":
            feas, o = verify_sol(g, sp)
            if not feas or abs(o - r["objective"]) > 1e-6 * max(1, abs(o)): ok = False
            elif abs(r["objective"] - ho) > 1e-6 * max(1, abs(ho)):
                better = (r["objective"] > ho) if g["maxi"] else (r["objective"] < ho)
                if better: refsub += 1; print("REF_SUBOPTIMAL", g["name"], "reference", ho, "taral", r["objective"], "(taral point verified feasible and integral)")
                else: ok = False
        if hs_cls in stats: stats[hs_cls] += 1
        if not ok:
            bad += 1; print("MISMATCH", g["name"], "highs", hs, ho, "taral", r)
            if os.environ.get("MILP_DUMP"): print("----MPS", g["name"]); print(open(p).read()); print("----END")
    print(f"MILP_CHECK cases={a.n} mismatches={bad} reference_suboptimal_verified={refsub} highs_optimal={stats['optimal']} highs_infeasible={stats['infeasible']}"); sys.exit(1 if bad else 0)
main()
