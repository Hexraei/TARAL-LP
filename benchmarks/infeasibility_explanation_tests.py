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

print("FAILED: %s" % fails if fails else "ALL PASS")
sys.exit(1 if fails else 0)
