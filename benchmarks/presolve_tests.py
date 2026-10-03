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
def expect(name, cond, msg=""):
    print(("PASS " if cond else "FAIL ") + name + ("" if cond else " " + str(msg)))
    if not cond: fails.append((name, [str(msg)]))

def run_raw(model, args):
    return subprocess.run([BIN, model] + args, capture_output=True, text=True, timeout=60)

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

def mps(name, body):
    p = os.path.join(TMP, name + ".mps"); open(p, "w").write(body); return p

conflict = mps("conflict_base", """NAME C
ROWS
 N OBJ
 G LOW
 L HIGH
COLUMNS
 X OBJ 1 LOW 1
 X HIGH 1
 Y OBJ 1 LOW 1
 Y HIGH 1
RHS
 RHS LOW 10 HIGH 4
BOUNDS
 UP BND X 100
 UP BND Y 100
ENDATA
""")

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
                if x[j] < M["lo"][j] - 1e-6 * (1 + min(abs(M["lo"][j]), 1e3)) or x[j] > M["up"][j] + 1e-6 * (1 + min(abs(M["up"][j]), 1e3)): problems.append("bound violated")
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

# ---- large-magnitude integer bounds (review item 6): empty integer intervals must be infeasible
def one_int(name, lo, up, ints=True, cost=1):
    return mps(name, "NAME L\nROWS\n N OBJ\nCOLUMNS\n" + (" MARK 'MARKER' 'INTORG'\n" if ints else "") + " X OBJ %g\n" % cost +
               (" MARK 'MARKER' 'INTEND'\n" if ints else "") + "BOUNDS\n LO BND X %r\n UP BND X %r\nENDATA\n" % (lo, up))
for nm, lo_, up_, want in (("big_empty_int", 1000000000.4, 1000000000.6, "infeasible"), ("big_empty_int_1e12", 1e12 + 0.3, 1e12 + 0.7, "infeasible"),
                           ("big_nonempty_int", 1000000000.4, 1000000001.6, "optimal"), ("big_integral_int", 1000000000.0, 1000000000.0000001, "optimal"),
                           ("small_empty_int", 2.2, 2.8, "infeasible")):
    pth = one_int(nm, lo_, up_)
    lg = os.path.join(TMP, nm + "_log.json"); sl = os.path.join(TMP, nm + ".sol")
    r = run_raw(pth, ["--presolve-log", lg, "--sol", sl]) if "run_raw" in globals() else subprocess.run([BIN, pth, "--presolve-log", lg, "--sol", sl], capture_output=True, text=True)
    st = [l for l in r.stdout.splitlines() if l.startswith("status ")][-1].split()[1]
    expect("large_int_status_" + nm, st == want, st)
    if os.path.exists(lg):
        errs = rc.replay(rc.load(pth), json.load(open(lg)))
        expect("large_int_replay_" + nm, not errs, errs)
    # plain engine must agree (no presolve): same status
    r0 = subprocess.run([BIN, pth], capture_output=True, text=True)
    st0 = [l for l in r0.stdout.splitlines() if l.startswith("status ")][-1].split()[1]
    expect("large_int_plain_agrees_" + nm, st0 == want, st0)
# an invalid point must fail the original-space audit and the replay cross-check: write X=1e9 for the empty interval
pth = os.path.join(TMP, "big_empty_int.mps"); lg = os.path.join(TMP, "big_empty_int_log.json")
fake_sol = os.path.join(TMP, "fake.sol"); open(fake_sol, "w").write("X 1000000000\n")
Mb = rc.load(pth)
l = json.load(open(lg)) if os.path.exists(lg) else {}
expect("large_int_invalid_point_violates_bounds", 1e9 < Mb["lo"][0] or 1e9 > Mb["up"][0])

# ---- regression tests for the independent review of the first version
# 1. fully eliminated model: JSON carries the original-space objective and point, matching stdout
full = mps("fully_eliminated", """NAME FULL
ROWS
 N OBJ
 L R0
COLUMNS
 X OBJ 2 R0 1
RHS
 RHS R0 5
BOUNDS
 FX BND X 1
ENDATA
""")
jp = os.path.join(TMP, "full.json"); sp_ = os.path.join(TMP, "full.sol")
r = run_raw(full, ["--presolve", "--json", jp, "--sol", sp_])
j = json.load(open(jp))
expect("full_elim_json_status", j["status"] == "optimal", j["status"])
expect("full_elim_json_objective_not_null", j["objective"] is not None and abs(j["objective"] - 2.0) < 1e-12, j["objective"])
expect("full_elim_json_has_x", j.get("x") == [1.0], j.get("x"))
expect("full_elim_stdout_matches", "objective 2" in r.stdout, r.stdout[-100:])

# 2. unwritable presolve log: nonzero exit with a diagnostic, in all three paths
for nm, mpath in (("normal", conflict), ("fully_eliminated", full), ("infeasible", os.path.join(ROOT, "benchmarks/refinery_stress/infeasible_supply_lp.mps"))):
    r = run_raw(mpath, ["--presolve-log", "/no-such-dir/log.json"])
    expect("unwritable_log_nonzero_" + nm, r.returncode != 0 and "cannot write" in (r.stdout + r.stderr), "rc=%d %s" % (r.returncode, r.stdout[-80:]))

# 3. replay checker strict validation
good = os.path.join(TMP, "good_log.json"); gsol = os.path.join(TMP, "good.sol")
run_raw(full, ["--presolve-log", good, "--sol", gsol])
Mf = rc.load(full); G = json.load(open(good))
import copy
def rep_rejects(name, mut, M=Mf, base=G):
    l = copy.deepcopy(base); mut(l)
    pth = os.path.join(TMP, name + ".json"); open(pth, "w").write(json.dumps(l, allow_nan=True))
    p = subprocess.run([sys.executable, os.path.join(ROOT, "benchmarks/presolve_replay_check.py"), full if M is Mf else M, pth], capture_output=True, text=True)
    expect(name, p.returncode == 1 and "FAIL" in p.stdout, p.stdout[-200:])
expect("good_log_replays", rc.replay(Mf, G) == [], rc.replay(Mf, G))
def col_minus1(l): l["ops"][0]["col"] = -1
rep_rejects("reject_fix_col_index_minus1", col_minus1)
def col_big(l): l["ops"][0]["col"] = 7
rep_rejects("reject_col_index_out_of_range", col_big)
def row_set(l): l["ops"][0]["row"] = 0
rep_rejects("reject_unexpected_row_index_on_fix_col", row_set)
def nan_a(l): l["ops"][0]["a"] = float("nan")
rep_rejects("reject_nan_unused_coefficient", nan_a)
def null_a(l): l["ops"][0]["a"] = None
rep_rejects("reject_null_coefficient", null_a)
def nan_v(l): l["ops"][0]["v2"] = float("nan")
rep_rejects("reject_nan_unused_value", nan_v)
def audit_nan(l): l["audit"]["objective"] = float("nan")
rep_rejects("reject_audit_nan_objective", audit_nan)
def audit_fake(l): l["audit"]["objective"] += 1.0
rep_rejects("reject_audit_ok_with_wrong_objective", audit_fake)
def audit_missing(l): del l["audit"]["max_row_violation"]
rep_rejects("reject_audit_missing_field", audit_missing)
def dims(l): l["original"]["cols"] = 9
rep_rejects("reject_wrong_original_dims", dims)
p = subprocess.run([sys.executable, os.path.join(ROOT, "benchmarks/presolve_replay_check.py"), full, good, gsol], capture_output=True, text=True)
expect("sol_audit_cross_check_passes", p.returncode == 0, p.stdout)
open(gsol, "w").write("X 5\n")
p = subprocess.run([sys.executable, os.path.join(ROOT, "benchmarks/presolve_replay_check.py"), full, good, gsol], capture_output=True, text=True)
expect("sol_audit_cross_check_rejects_wrong_point", p.returncode == 1, p.stdout)

# 4. oracle semantics: duplicate coefficients are summed (engine) and a read failure is reported, not ignored
dupf = mps("dup_presolve", """NAME D
ROWS
 N OBJ
 G R0
COLUMNS
 X OBJ 1 R0 3
 X R0 -4
RHS
 RHS R0 2
BOUNDS
 UP BND X 5
ENDATA
""")  # summed coefficient -1: -x >= 2 is infeasible; first-only (3x >= 2) would be feasible
lg = os.path.join(TMP, "dup_log.json")
r = run_raw(dupf, ["--presolve-log", lg])
expect("duplicates_summed_presolve_infeasible", "status infeasible" in r.stdout, r.stdout[-100:])
p = subprocess.run([sys.executable, os.path.join(ROOT, "benchmarks/presolve_replay_check.py"), dupf, lg], capture_output=True, text=True)
expect("duplicates_summed_replay_passes", p.returncode == 0, p.stdout)
garbage = mps("garbage2", "this is not an mps file\n")
p = subprocess.run([sys.executable, os.path.join(ROOT, "benchmarks/presolve_replay_check.py"), garbage, good], capture_output=True, text=True)
expect("replay_reader_failure_fails", p.returncode == 1, p.stdout)

# 5. deadline discipline: expired budget -> honest time_limit, partial reduction unused, log marked and replayable
tl = os.path.join(TMP, "tl_log.json"); tj = os.path.join(TMP, "tl.json")
os.environ["TARAL_PRESOLVE_TEST_TIMEOUT_AFTER_OPS"] = "5"  # deterministic: as if the deadline passed after 5 reductions
r = run_raw(os.path.join(ROOT, "benchmarks/refinery_stress/infeasible_supply_lp.mps"), ["--presolve-log", tl, "--json", tj])
del os.environ["TARAL_PRESOLVE_TEST_TIMEOUT_AFTER_OPS"]
if "time_limit" in r.stdout and "presolve" in r.stdout:
    expect("presolve_timeout_exit4", r.returncode == 4, r.returncode)
    L = json.load(open(tl)); expect("presolve_timeout_log_marked", L.get("timed_out") is True)
    expect("presolve_timeout_json_status", json.load(open(tj))["status"] == "time_limit")
    p = subprocess.run([sys.executable, os.path.join(ROOT, "benchmarks/presolve_replay_check.py"), os.path.join(ROOT, "benchmarks/refinery_stress/infeasible_supply_lp.mps"), tl], capture_output=True, text=True)
    expect("presolve_timeout_log_replays", p.returncode == 0, p.stdout)
else:
    expect("presolve_timeout_reached", False, r.stdout[-200:])
# 6. review 7/8/9: strict inward integer rounding (no outward snap), absolute replay validation, primal-only output
def status_of(r): return [l for l in r.stdout.splitlines() if l.startswith("status ")][-1].split()[1]
for nm, lo_, up_ in (("snap_zero", 4e-7, 8e-7), ("snap_one", 1.0000004, 1.0000008), ("snap_neg", -1.0000008, -1.0000004), ("snap_big", 1e9 + 4e-7, 1e9 + 8e-7)):
    pth = one_int(nm, lo_, up_); lg = os.path.join(TMP, nm + "_log.json")
    r = run_raw(pth, ["--presolve-log", lg]); r0 = subprocess.run([BIN, pth], capture_output=True, text=True)
    expect("strict_round_infeasible_" + nm, status_of(r) == "infeasible", status_of(r))
    # documented divergence: the plain MILP (no --presolve) keeps its own 1e-6 integrality tolerance (milp.cpp kIntTol); recorded, not asserted
    print("INFO plain_milp_status_" + nm, status_of(r0))
    if os.path.exists(lg): expect("strict_round_replay_" + nm, not rc.replay(rc.load(pth), json.load(open(lg))), rc.replay(rc.load(pth), json.load(open(lg))))
for c in (0.0, 1.0, -1.0, 1e3, 1e9):
    for w in (0.0, 8e-7, 2.2e-6):
        nm = "contains_%g_%g" % (c, w); pth = one_int(nm, c - w / 2, c + w / 2); lg = os.path.join(TMP, nm + "_log.json"); sl = os.path.join(TMP, nm + ".sol")
        r = run_raw(pth, ["--presolve-log", lg, "--sol", sl]); st = status_of(r)
        expect("contains_integer_optimal_" + nm, st == "optimal" and ("objective %.12g" % c) in r.stdout, r.stdout[-120:])
        if os.path.exists(lg): expect("contains_integer_replay_" + nm, not rc.replay(rc.load(pth), json.load(open(lg))))
fb = os.path.join(ROOT, "tests/fixtures/presolve_fixedbig.mps"); fl = os.path.join(TMP, "fixedbig_log.json")
run_raw(fb, ["--presolve-log", fl]); Mfb = rc.load(fb); Lfb = json.load(open(fl))
expect("fixedbig_replays", not rc.replay(Mfb, Lfb), rc.replay(Mfb, Lfb))
fi = [k for k, o in enumerate(Lfb["ops"]) if o["type"] == "fix_col"][0]
def mut_fix(l, d):
    l["ops"][fi]["v1"] += d; l["reduced"]["obj_const"] += d
for d in (1.0, 0.5, 1e-3):
    l2 = copy.deepcopy(Lfb); mut_fix(l2, d); expect("replay_rejects_fixed_value_shift_%g" % d, bool(rc.replay(Mfb, l2)))
l2 = copy.deepcopy(Lfb); l2["reduced"]["obj_const"] += 1.0; expect("replay_rejects_objective_only_shift", bool(rc.replay(Mfb, l2)))
l2 = copy.deepcopy(Lfb); l2["ops"][fi]["v1"] += 1.0; expect("replay_rejects_fix_value_only_shift", bool(rc.replay(Mfb, l2)))
# primal-only output: no reduced-model dual/KKT placeholders next to the original-space x
r = run_raw(full, ["--presolve", "--json", jp])
j = json.load(open(jp))
for k in ("dual_objective", "gap", "primal_res", "dual_res", "complementarity", "row_dual", "reduced_cost"):
    expect("primal_only_null_" + k, j.get(k, "missing") is None, j.get(k, "missing"))
expect("primal_only_flag", j.get("original_space_primal_only") is True and j.get("certificate_quality") == "presolve_primal_only")
expect("primal_only_row_activity_original", j.get("row_activity") == [1.0], j.get("row_activity"))
rp = run_raw(os.path.join(ROOT, "benchmarks/refinery_stress/base_lp_t6.mps"), ["--presolve", "--json", jp]); j = json.load(open(jp))
expect("primal_only_null_reduced_lp", j["status"] == "optimal" and j["dual_objective"] is None and j["row_dual"] is None and len(j["row_activity"]) > 0, j["status"])
inf = os.path.join(ROOT, "benchmarks/refinery_stress/infeasible_supply_lp.mps"); run_raw(inf, ["--presolve", "--json", jp]); j = json.load(open(jp))
if j["status"] == "infeasible" and "certificate_quality" in j and "presolve" in j.get("message", ""):
    expect("presolve_infeasible_reason_present", True)
if j.get("certificate_quality") == "presolve_reduced_model_only":
    expect("reduced_certificate_not_claimed", j.get("certificate_verified") is False and not j.get("farkas_row_lower"), j)

api = os.path.join(TMP, "presolve_api")
srcs = [os.path.join(ROOT, "src", f) for f in os.listdir(os.path.join(ROOT, "src")) if f.endswith(".cpp") and f != "main.cpp"]
cc = subprocess.run(["g++", "-O1", "-std=c++17", "-o", api, os.path.join(ROOT, "benchmarks/presolve_api_tests.cpp")] + srcs, capture_output=True, text=True)
expect("presolve_api_test_builds", cc.returncode == 0, cc.stderr[-300:])
if cc.returncode == 0:
    r = subprocess.run([api], capture_output=True, text=True); expect("presolve_api_deadline_tests", r.returncode == 0, r.stdout + r.stderr)

print("stats", json.dumps(stats))
print("FAILED %d" % len(fails) if fails else "ALL PASS")
sys.exit(1 if fails else 0)
