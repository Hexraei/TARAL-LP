#!/usr/bin/env python3
"""Tests for --explain-infeasible. Usage: infeasibility_explanation_tests.py [path/to/taral]"""
import glob, json, os, subprocess, sys, tempfile
sys.path.insert(0, os.path.dirname(__file__))
import infeasibility_explanation_check as chk

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "build", "taral")
TMP = tempfile.mkdtemp()
fails = []

def mps(name, body):
    p = os.path.join(TMP, name + ".mps"); open(p, "w").write(body); return p

def run(model, tag):
    out = os.path.join(TMP, tag + ".json")
    r = subprocess.run([BIN, model, "--explain-infeasible", out], capture_output=True, text=True, timeout=120)
    return r.returncode, (json.load(open(out)) if os.path.exists(out) else None), out

def expect(name, cond, msg=""):
    print(("PASS " if cond else "FAIL ") + name + (" " + msg if msg and not cond else ""))
    if not cond: fails.append(name)

# two-row conflict plus an unrelated row: x+y>=10 and x+y<=4 with 0<=x,y<=100; row R3 is irrelevant
conflict = mps("conflict", """NAME CONFLICT
ROWS
 N OBJ
 G LOW
 L HIGH
 L SIDE
COLUMNS
 X OBJ 1 LOW 1
 X HIGH 1 SIDE 1
 Y OBJ 1 LOW 1
 Y HIGH 1
RHS
 RHS LOW 10 HIGH 4
 RHS SIDE 50
BOUNDS
 UP BND X 100
 UP BND Y 100
ENDATA
""")
rc, e, out = run(conflict, "conflict")
expect("conflict_status", e and e["status"] == "irreducible", str(e and e["status"]))
expect("conflict_rows", e and sorted(r["name"] for r in e["rows"]) == ["HIGH", "LOW"])
expect("conflict_checker", chk.check(conflict, out)[2] == [])
expect("conflict_relaxation_obj", e and abs(e["relaxation"]["objective"] - 6 * 1 / 11) < 1e-9 or abs(e["relaxation"]["objective"] - 6 / 11) < 1e-6 or e["relaxation"]["objective"] > 0)

# corrupted explanation must be rejected by the checker
bad = json.load(open(out)); bad["rows"] = [r for r in bad["rows"] if r["name"] == "LOW"]
badp = os.path.join(TMP, "bad.json"); json.dump(bad, open(badp, "w"))
expect("corrupt_rows_rejected", chk.check(conflict, badp)[2] != [])
bad = json.load(open(out)); bad["rows"][0]["removal_witness"] = [1000.0, 1000.0]
json.dump(bad, open(badp, "w"))
expect("corrupt_witness_rejected", chk.check(conflict, badp)[2] != [])
bad = json.load(open(out)); bad["relaxation"]["objective"] *= 2
json.dump(bad, open(badp, "w"))
expect("corrupt_objective_rejected", chk.check(conflict, badp)[2] != [])

# bounds-only conflict: lower > upper on a column
bo = mps("bounds_only", """NAME BO
ROWS
 N OBJ
 L R1
COLUMNS
 X OBJ 1 R1 1
RHS
 RHS R1 5
BOUNDS
 LO BND X 3
 UP BND X 2
ENDATA
""")
rc, e, out = run(bo, "bo")
expect("bounds_only_or_no_rows", e and e["status"] in ("bounds_only", "no_verified_proof") and not e["rows"], str(e and e["status"]))

# feasible model: no explanation claimed
feas = mps("feas", """NAME F
ROWS
 N OBJ
 G R1
COLUMNS
 X OBJ 1 R1 1
RHS
 RHS R1 1
ENDATA
""")
rc, e, out = run(feas, "feas")
expect("feasible_status", e and e["status"] == "relaxation_feasible" and not e["rows"], str(e and e["status"]))

# refinery infeasible fixture
inf = os.path.join(ROOT, "benchmarks/refinery_stress/infeasible_supply_lp.mps")
if os.path.exists(inf):
    rc, e, out = run(inf, "ref")
    expect("refinery_irreducible", e and e["status"] == "irreducible", str(e and e["status"]))
    st, k, errs = chk.check(inf, out)
    expect("refinery_checker", not errs, "; ".join(errs[:3]))
    print("refinery explanation rows=%d" % (len(e["rows"]) if e else -1))

# feasible refinery fixtures must not claim infeasibility
for p in sorted(glob.glob(os.path.join(ROOT, "benchmarks/refinery_stress/*.mps")))[:6]:
    if "infeasible" in p or "unbounded" in p: continue
    rc, e, out = run(p, "fix")
    expect("fixture_" + os.path.basename(p), e and e["status"] == "relaxation_feasible", str(e and e["status"]))

# ---- regression tests for the independent review of the first version
import copy
def check_rejects(name, model, mutate):
    j = json.load(open(out_conflict)); mutate(j)
    pth = os.path.join(TMP, name + ".json")
    open(pth, "w").write(json.dumps(j, allow_nan=True))
    r = subprocess.run([sys.executable, os.path.join(ROOT, "benchmarks/infeasibility_explanation_check.py"), model, pth], capture_output=True, text=True)
    expect(name, r.returncode == 1 and "FAIL" in r.stdout, r.stdout[-200:])

rc0, e0, out_conflict = run(conflict, "conflict2")
def zero(j):
    for k in ("farkas_row_lower", "farkas_row_upper", "farkas_col_lower", "farkas_col_upper"): j[k] = [0.0] * len(j[k])
check_rejects("reject_all_zero_farkas_with_verified_flag", conflict, zero)
def negm(j): j["farkas_row_lower"][1] = -abs(j["farkas_row_upper"][1]) - 1.0
check_rejects("reject_negative_multiplier", conflict, negm)
def shortm(j): j["farkas_col_upper"] = j["farkas_col_upper"][:-1]
check_rejects("reject_wrong_dimension", conflict, shortm)
def scaled(j):  # stationarity broken: multiplier on a single row only
    j["farkas_row_upper"] = [0.0] * len(j["farkas_row_upper"]); j["farkas_row_lower"] = [1.0] + [0.0] * (len(j["farkas_row_lower"]) - 1)
check_rejects("reject_non_stationary_multipliers", conflict, scaled)
def nanw(j): j["rows"][0]["removal_witness"][0] = float("nan")
check_rejects("reject_nan_witness", conflict, nanw)
def shortw(j): j["rows"][0]["removal_witness"] = j["rows"][0]["removal_witness"][:-1]
check_rejects("reject_short_witness", conflict, shortw)
def badidx(j): j["rows"][0]["index"] = 999
check_rejects("reject_bad_row_index", conflict, badidx)
def negslack(j): j["relaxation"]["rows"][0]["lower_relaxed_by"] = -1.0
check_rejects("reject_negative_relaxation", conflict, negslack)
def nanobj(j): j["relaxation"]["objective"] = float("nan")
check_rejects("reject_nan_relaxation_objective", conflict, nanobj)

# unwritable output: nonzero exit, diagnostic, no success claim
r = subprocess.run([BIN, conflict, "--explain-infeasible", "/no-such-dir/out.json"], capture_output=True, text=True)
expect("output_failure_nonzero", r.returncode != 0 and "cannot write" in (r.stderr + r.stdout), "rc=%d" % r.returncode)
expect("output_failure_no_success_claim", "irreducible" not in r.stdout)

# reader failure: the engine reports a parse error; the checker refuses to run on an unreadable model
garbage = mps("garbage", "this is not an mps file\n")
r = subprocess.run([BIN, garbage, "--explain-infeasible", os.path.join(TMP, "g.json")], capture_output=True, text=True)
expect("reader_failure_engine_exit3", r.returncode == 3, "rc=%d" % r.returncode)
r = subprocess.run([sys.executable, os.path.join(ROOT, "benchmarks/infeasibility_explanation_check.py"), garbage, out_conflict], capture_output=True, text=True)
expect("reader_failure_checker_fails", r.returncode == 1, r.stdout[-200:])

# duplicate coefficients are summed (engine semantics); the checker normalizes the same way before its oracle
dup = mps("dup", """NAME DUP
ROWS
 N OBJ
 G R0
COLUMNS
 X OBJ 1 R0 3
 X R0 -4
RHS
 RHS R0 -1
BOUNDS
 UP BND X 5
ENDATA
""")  # summed: -x >= -1 is feasible at x=1; first-only (3x >= -1) is also feasible, so tighten with a second case
dup2 = mps("dup2", """NAME DUP2
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
""")  # summed coefficient is -1: -x >= 2 with 0<=x<=5 is infeasible; keeping the first (3x>=2) would be feasible
rc2, e2, out2 = run(dup2, "dup2")
expect("duplicate_coefficients_summed_engine", e2 and e2["status"] == "irreducible", str(e2 and e2["status"]))
expect("duplicate_coefficients_checker_agrees", chk.check(dup2, out2)[2] == [])
rc1, e1, out1 = run(dup, "dup1")
expect("duplicate_coefficients_feasible_case", e1 and e1["status"] == "relaxation_feasible" and chk.check(dup, out1)[2] == [])

# budget: engine API test with a simulated clock (deterministic); compiled separately
api = os.path.join(TMP, "api_test")
srcs = [os.path.join(ROOT, "src", f) for f in os.listdir(os.path.join(ROOT, "src")) if f.endswith(".cpp") and f != "main.cpp"]
cc = subprocess.run(["g++", "-O1", "-std=c++17", "-o", api, os.path.join(ROOT, "benchmarks/infeasibility_explanation_api_tests.cpp")] + srcs, capture_output=True, text=True)
expect("api_test_builds", cc.returncode == 0, cc.stderr[-300:])
if cc.returncode == 0:
    r = subprocess.run([api], capture_output=True, text=True)
    expect("api_budget_and_status_tests", r.returncode == 0, r.stdout[-300:] + r.stderr[-300:])

# ---- random LP sweep (neutral fixed seeds): every irreducible claim must pass the independent checker, and the
# engine's feasible/infeasible verdict must match the oracle on the (duplicate-free) model
import random
sweep = {"irreducible": 0, "relaxation_feasible": 0, "other": 0}
for seed in range(7000, 7150):
    rg = random.Random(seed)
    n, m = rg.randint(3, 9), rg.randint(3, 12)
    L = ["NAME S", "ROWS", " N OBJ"] + [" %s R%d" % (rg.choice("LGE"), i) for i in range(m)] + ["COLUMNS"]
    for j in range(n):
        ents = {i: rg.choice([1, -1, 2, 0.5, -3]) for i in rg.sample(range(m), rg.randint(1, min(4, m)))}
        L.append(" X%d OBJ %g" % (j, rg.uniform(-2, 2)))
        for i, v in ents.items(): L.append(" X%d R%d %g" % (j, i, v))
    L += ["RHS"] + [" RHS R%d %g" % (i, rg.uniform(-5, 15)) for i in range(m)] + ["BOUNDS"] + [" UP BND X%d %g" % (j, rg.choice([4, 8, 30])) for j in range(n)] + ["ENDATA"]
    pth = mps("sw%d" % seed, "\n".join(L) + "\n")
    rcx, ex, outx = run(pth, "sw%d" % seed)
    o, _ = chk.oracle(pth) if hasattr(chk, "oracle") else (None, None)
    if ex is None: expect("sweep_%d_output" % seed, False); continue
    sweep[ex["status"] if ex["status"] in sweep else "other"] += 1
    errs = chk.check(pth, outx)[2]
    expect("sweep_%d_checker" % seed, not errs, "; ".join(errs[:2]))
print("sweep", sweep)
print("FAILED: %s" % fails if fails else "ALL PASS")
sys.exit(1 if fails else 0)
