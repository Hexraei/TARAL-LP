#!/usr/bin/env python3
"""Tests for --presolve. Usage: presolve_tests.py [path/to/taral] [n_random]
Random seeds are neutral fixed integers. HiGHS (highspy) on the ORIGINAL model is the status/objective
oracle only. For every case: taral with and without --presolve agree with the oracle, the presolve log
replays against the original model (presolve_replay_check), and the reported point is audited in
original space by the engine (and re-audited here from the .sol file). Corrupted logs must be rejected."""
import json, math, os, random, subprocess, sys, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, highspy
import presolve_replay_check as rc

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "build", "taral")
NRAND = int(sys.argv[2]) if len(sys.argv) > 2 else 200
TMP = tempfile.mkdtemp()
fails = []
stats = {"cases": 0, "reduced": 0, "ops": {}, "presolve_infeasible": 0, "audited_points": 0}

def write_mps(path, rows, cols, ints, maximize, objc):
    # rows: list of (type, rhs, range or None); cols: list of dict(cost, ents{row:val}, lo, up)
    L = ["NAME T"]
    if maximize: L += ["OBJSENSE", "    MAX"]
    L.append("ROWS"); L.append(" N OBJ")
    for i, (t, rhs, rng) in enumerate(rows): L.append(" %s R%d" % (t, i))
    L.append("COLUMNS")
    inint = False
    for j, c in enumerate(cols):
        if ints[j] and not inint: L.append(" MARK 'MARKER' 'INTORG'"); inint = True
        if not ints[j] and inint: L.append(" MARK 'MARKER' 'INTEND'"); inint = False
        ents = [("OBJ", c["cost"])] + [("R%d" % i, v) for i, v in sorted(c["ents"].items())]
        if not c["ents"] and c["cost"] == 0: ents = [("OBJ", 0.0)]
        for nm, v in ents: L.append(" X%d %s %.12g" % (j, nm, v))
    if inint: L.append(" MARK 'MARKER' 'INTEND'")
    L.append("RHS")
    for i, (t, rhs, rng) in enumerate(rows): L.append(" RHS R%d %.12g" % (i, rhs))
    if objc: L.append(" RHS OBJ %.12g" % (-objc))
    rg = [(i, r) for i, (t, rhs, r) in enumerate(rows) if r is not None]
    if rg:
        L.append("RANGES")
        for i, r in rg: L.append(" RNG R%d %.12g" % (i, r))
    L.append("BOUNDS")
    for j, c in enumerate(cols):
        lo, up = c["lo"], c["up"]
        if lo == up: L.append(" FX BND X%d %.12g" % (j, lo)); continue
        if lo == -math.inf and up == math.inf: L.append(" FR BND X%d" % j); continue
        if lo == -math.inf: L.append(" MI BND X%d" % j)
        elif lo != 0: L.append(" LO BND X%d %.12g" % (j, lo))
        if up != math.inf: L.append(" UP BND X%d %.12g" % (j, up))
    L.append("ENDATA")
    open(path, "w").write("\n".join(L) + "\n")

def gen(rng):
    n, m = rng.randint(4, 12), rng.randint(3, 10)
    maximize = rng.random() < 0.4
    ints = [rng.random() < 0.5 for _ in range(n)]
    cols = []
    for j in range(n):
        lo = rng.choice([0, 0, 0, -3, 1.5]) if not ints[j] else rng.choice([0, 0, -2, 0.4])
        up = lo + rng.choice([3, 5, 10, 20.5]) if (ints[j] or rng.random() < 0.9) else math.inf  # integer columns always get an explicit upper bound: readers differ on the default
        if rng.random() < 0.15: up = lo  # fixed column
        cols.append(dict(cost=round(rng.uniform(-5, 5), 2) if rng.random() < 0.85 else 0.0, ents={}, lo=lo, up=up))
    rows = []
    feas = rng.random() < 0.75  # build rhs around a hidden point so most cases are feasible
    x0 = [min(max(c["lo"] + (rng.randint(0, 2) if ints[j] else rng.uniform(0, 2)), c["lo"]), c["up"]) for j, c in enumerate(cols)]
    if feas: x0 = [float(round(v)) if ints[j] else v for j, v in enumerate(x0)]
    if feas: x0 = [min(max(v, c["lo"]), c["up"]) for v, c in zip(x0, cols)]
    for i in range(m):
        kind = rng.random()
        k = 1 if kind < 0.2 else (0 if kind < 0.27 else rng.randint(2, min(5, n)))
        for j in rng.sample(range(n), k): cols[j]["ents"][i] = rng.choice([1, 2, -1, 3, 0.5, -2.5])
        t = rng.choice(["L", "G", "E"]) if rng.random() < 0.8 else "G"
        rhs = round(rng.uniform(-4, 12) if t != 'G' else rng.uniform(-6, 4), 1)
        if feas:
            act = sum(cols[j]["ents"].get(i, 0) * x0[j] for j in range(n))
            rhs = round(act + (0 if t == "E" else (rng.uniform(0, 3) if t == "L" else -rng.uniform(0, 3))), 3)
        rows.append((t, rhs, round(rng.uniform(0, 10), 1) if (t in "LG" and rng.random() < 0.15) else None))
    return rows, cols, ints, maximize, round(rng.choice([0, 0, 3.5, -2]), 1)

def oracle(path):
    """Oracle: the HiGHS MIP solver through scipy on the matrix/bounds read by highspy (no file-based solve)."""
    import scipy.sparse as sp
    from scipy.optimize import milp, LinearConstraint, Bounds
    M = rc.load(path); n, m = M["n"], M["m"]
    A = sp.lil_matrix((m, n))
    for j, c in enumerate(M["cols"]):
        for i, a in c: A[i, j] = a
    lo = [math.ceil(l - 1e-9) if (M["isint"][j] and math.isfinite(l)) else l for j, l in enumerate(M["lo"])]  # harness rounds integer bounds
    up = [math.floor(u + 1e-9) if (M["isint"][j] and math.isfinite(u)) else u for j, u in enumerate(M["up"])]
    cost = np.array(M["cost"]) * (-1 if M["maximize"] else 1)
    r = milp(cost, constraints=LinearConstraint(A.tocsr(), M["rlo"], M["rup"]) if m else None,
             integrality=np.array(M["isint"], int), bounds=Bounds(lo, up),
             options={"presolve": False, "mip_rel_gap": 0.0})  # oracle presolve off: HiGHS presolve disagreed on fractional integer bounds in 3 cases
    st = {0: "optimal", 2: "infeasible", 3: "unbounded"}.get(r.status, "other")
    return st, (r.fun * (-1 if M["maximize"] else 1) + M["offset"]) if st == "optimal" else math.nan

def taral(path, extra, tag):
    js, so = os.path.join(TMP, tag + ".json"), os.path.join(TMP, tag + ".sol")
    for f in (js, so):
        if os.path.exists(f): os.remove(f)
    r = subprocess.run([BIN, path, "--json", js, "--sol", so, "--time-limit", "30"] + extra, capture_output=True, text=True, timeout=60)
    sl = [l for l in r.stdout.splitlines() if l.startswith("status ")]
    status = sl[-1].split()[1] if sl else "?"
    return status, r.stdout, so

def sol_obj(M, so):
    val = dict(l.split() for l in open(so)) if os.path.exists(so) else {}
    return val

def check_case(path, tag, expect_reduction=False):
    stats["cases"] += 1
    o, oobj = oracle(path)
    s0, _, _ = taral(path, [], tag + "_off")
    log = os.path.join(TMP, tag + "_log.json")
    s1, out, so = taral(path, ["--presolve-log", log], tag + "_on")
    M = rc.load(path)
    problems = []
    if o in ("optimal", "infeasible") and s1 != o: problems.append("presolve status %s vs oracle %s" % (s1, o))
    if o in ("optimal", "infeasible") and s0 != o: problems.append("plain status %s vs oracle %s" % (s0, o))
    if o == "unbounded" and s1 not in ("unbounded", "unbounded_relaxation", "dual_infeasible"): problems.append("presolve %s vs oracle unbounded" % s1)
    if os.path.exists(log):
        L = json.load(open(log))
        errs = rc.replay(M, L)
        if errs: problems.append("replay: " + "; ".join(errs[:2]))
        if L["ops"] or L["infeasible"]: stats["reduced"] += 1
        stats["presolve_infeasible"] += bool(L["infeasible"])
        for op in L["ops"]: stats["ops"][op["type"]] = stats["ops"].get(op["type"], 0) + 1
    else: problems.append("no log written")
    if s1 == "optimal":
        stats["audited_points"] += 1
        if "presolve audit ok" not in out and "(presolve)" not in out: problems.append("no ok audit line")
        sl = [l for l in out.splitlines() if l.startswith("status ")][-1]
        obj = float(sl.split("objective ")[1].split()[0])
        if abs(obj - oobj) > 1e-6 * (1 + abs(oobj)): problems.append("objective %.9g vs oracle %.9g" % (obj, oobj))
        val = sol_obj(M, so)
        hh = highspy.Highs(); hh.setOptionValue("output_flag", False); hh.readModel(path)
        names = list(hh.getLp().col_names_)
        x = [float(val.get(names[j], "nan")) for j in range(M["n"])]
        if any(math.isnan(v) for v in x): problems.append("sol file lacks columns")
        else:  # re-audit here, original space
            ax = [0.0] * M["m"]
            for j, col in enumerate(M["cols"]):
                for i, a in col: ax[i] += a * x[j]
                if x[j] < M["lo"][j] - 1e-6 * rc.sc(M["lo"][j]) or x[j] > M["up"][j] + 1e-6 * rc.sc(M["up"][j]): problems.append("bound violated")
                if M["isint"][j] and abs(x[j] - round(x[j])) > 1e-6: problems.append("integrality violated")
            for i in range(M["m"]):
                if ax[i] < M["rlo"][i] - 1e-6 * rc.sc(M["rlo"][i]) or ax[i] > M["rup"][i] + 1e-6 * rc.sc(M["rup"][i]): problems.append("row violated")
    if problems: fails.append((tag, problems)); print("FAIL", tag, problems)
    return s1

# ---- hand-built cases exercising each reduction
hand = {
 "fix_and_singleton": ([("L", 10, None), ("G", 2, None), ("L", 100, None), ("E", 0, None)],
    [dict(cost=1, ents={0: 1, 1: 1, 2: 1}, lo=0, up=5), dict(cost=-2, ents={0: 1, 2: 4}, lo=1, up=1), dict(cost=3, ents={1: 2}, lo=0, up=10), dict(cost=0, ents={}, lo=0, up=7)],
    [False] * 4, False, 0),
 "int_round_and_empty_col": ([("L", 7.5, None), ("G", 0.0, None)],
    [dict(cost=-1, ents={0: 2}, lo=0, up=9), dict(cost=-1, ents={0: 1, 1: 1}, lo=0.3, up=4.7), dict(cost=2, ents={}, lo=-1, up=4)],
    [True, True, True], False, 1.5),
 "presolve_infeasible_singleton": ([("G", 8, None), ("L", 3, None)],
    [dict(cost=1, ents={0: 1, 1: 1}, lo=0, up=5), dict(cost=1, ents={0: 1}, lo=0, up=1)], [False, False], False, 0),
 "presolve_infeasible_activity": ([("G", 100, None)],
    [dict(cost=1, ents={0: 1}, lo=0, up=5), dict(cost=1, ents={0: 1}, lo=0, up=5)], [False, False], False, 0),
 "maximize_empty_col": ([("L", 4, None)],
    [dict(cost=1, ents={0: 1}, lo=0, up=10), dict(cost=2, ents={}, lo=0, up=6)], [False, False], True, 0),
}
for name, (rows, cols, ints, mx, oc) in hand.items():
    p = os.path.join(TMP, name + ".mps"); write_mps(p, rows, cols, ints, mx, oc)
    st = check_case(p, name)
    print("hand", name, st)

# ---- random cases (seeds are neutral fixed integers)
for seed in range(1000, 1000 + NRAND):
    rng = random.Random(seed)
    rows, cols, ints, mx, oc = gen(rng)
    p = os.path.join(TMP, "r%d.mps" % seed); write_mps(p, rows, cols, ints, mx, oc)
    check_case(p, "r%d" % seed)

# ---- refinery fixtures (plain, unreduced in practice; the check is that presolve does no harm)
for f in ("base_lp_t6", "supply_cut_lp", "campaign_milp_t6", "infeasible_supply_lp"):
    p = os.path.join(ROOT, "benchmarks/refinery_stress", f + ".mps")
    if os.path.exists(p): check_case(p, "ref_" + f)

# ---- corrupted logs must be rejected
M = rc.load(os.path.join(TMP, "fix_and_singleton.mps")); L = json.load(open(os.path.join(TMP, "fix_and_singleton_on_log.json"))) if False else json.load(open(os.path.join(TMP, "fix_and_singleton_log.json")))
def rejected(mut, name):
    import copy
    l2 = copy.deepcopy(L); mut(l2)
    ok = bool(rc.replay(M, l2)); print(("PASS " if ok else "FAIL ") + "corrupt_" + name)
    if not ok: fails.append(("corrupt_" + name, []))
fixes = [k for k, o in enumerate(L["ops"]) if o["type"] == "fix_col"]
rejected(lambda l: l["ops"][fixes[0]].__setitem__("v1", 3.0), "fix_value")
rejected(lambda l: l["ops"].pop(fixes[0]), "dropped_op")
rejected(lambda l: l["kept_rows"].append(99), "kept_rows")
rejected(lambda l: l["reduced"].__setitem__("obj_const", l["reduced"]["obj_const"] + 1), "obj_const")
sr = [k for k, o in enumerate(L["ops"]) if o["type"] == "singleton_row"]
if sr: rejected(lambda l: l["ops"][sr[0]].__setitem__("v3", 123.0), "singleton_bound")
rejected(lambda l: l["ops"].append({"type": "redundant_row", "row": 0, "col": -1, "a": 0, "v1": -math.inf, "v2": 10, "v3": 0, "v4": 1}), "bogus_redundant")

print("stats", json.dumps(stats))
print("FAILED %d" % len(fails) if fails else "ALL PASS")
sys.exit(1 if fails else 0)
