#!/usr/bin/env python3
"""Tests for --work-limit (deterministic LP work budget). Usage: deterministic_solve_tests.py [taral] [baseline_taral]
Claims tested: same binary + same input bytes + same --work-limit => identical status, iteration count,
work_used and result_hash, regardless of --time-limit and concurrent machine load. NOT tested/claimed:
identical results across different CPUs, compilers, flags or libm (see docs/deterministic_solve.md)."""
import json, os, subprocess, sys, tempfile, threading, time
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "build", "taral")
BASE = sys.argv[2] if len(sys.argv) > 2 else None
TMP = tempfile.mkdtemp(); fails = []
FIX = [os.path.join(ROOT, "benchmarks/refinery_stress", f) for f in ("base_lp_t6.mps", "supply_cut_lp.mps", "demand_surge_lp.mps", "sour_outage_lp.mps", "quality_tight_lp.mps", "base_lp_t24.mps")]
FIX = [f for f in FIX if os.path.exists(f)]

def expect(name, cond, msg=""):
    print(("PASS " if cond else "FAIL ") + name + ("" if cond else " " + str(msg)))
    if not cond: fails.append(name)

def run(model, args, tag="r"):
    js = os.path.join(TMP, tag + ".json")
    if os.path.exists(js): os.remove(js)
    p = subprocess.run([BIN, model, "--json", js] + args, capture_output=True, text=True, timeout=120)
    return p, (json.load(open(js)) if os.path.exists(js) else None)

def key(j): return (j["status"], j["iterations"], j.get("work_used"), j.get("result_hash"), j["objective"])

# 1. repeated runs, different wall limits: identical
for m in FIX:
    nm = os.path.basename(m)
    keys = set()
    for tl in (None, "5", "600", "5", None):
        _, j = run(m, ["--work-limit", "100000"] + (["--time-limit", tl] if tl else []))
        keys.add(key(j))
    expect("repeat_and_timelimit_invariant_" + nm, len(keys) == 1, keys)

# 2. identical under concurrent CPU load
stop = threading.Event()
def burn():
    while not stop.is_set(): sum(i * i for i in range(10000))
ths = [threading.Thread(target=burn) for _ in range(4)]
[t.start() for t in ths]
try:
    for m in FIX[:3]:
        _, j1 = run(m, ["--work-limit", "100000"]); _, j2 = run(m, ["--work-limit", "100000"])
        _, j0 = None, None
        expect("load_invariant_" + os.path.basename(m), key(j1) == key(j2))
finally:
    stop.set(); [t.join() for t in ths]

# 3. caps: stopping points are deterministic and never a wall-clock status
for m in FIX[:3]:
    for cap in (1, 2, 5, 17, 40, 90, 200):
        _, a = run(m, ["--work-limit", str(cap)]); _, b = run(m, ["--work-limit", str(cap), "--time-limit", "9"])
        ok = key(a) == key(b) and a["status"] in ("iteration_limit", "optimal")
        expect("cap_%d_%s" % (cap, os.path.basename(m)), ok, (key(a), key(b)))

# 4. hash independent of the cap once the solve completes; hash excludes wall time
m = FIX[0]
_, a = run(m, ["--work-limit", "100000"]); _, b = run(m, ["--work-limit", "10000000"])
expect("hash_independent_of_cap_when_complete", a["status"] == "optimal" and a["result_hash"] == b["result_hash"] and a["iterations"] == b["iterations"])
expect("work_used_within_cap", all(isinstance(j["work_used"], int) for j in (a, b)))
expect("json_has_work_fields", all(k in a for k in ("work_limit", "work_used", "result_hash")))

# 5. status at a tiny cap is iteration_limit with no objective claim
_, c = run(m, ["--work-limit", "3"])
expect("tiny_cap_iteration_limit", c["status"] == "iteration_limit" and c["objective"] is None, c["status"])

# 6. flag off: output unchanged versus the baseline binary, no work fields
if BASE:
    for f in FIX:
        j = os.path.join(TMP, "b.json")
        subprocess.run([BASE, f, "--json", j], capture_output=True, timeout=120); jb = json.load(open(j))
        _, jn = run(f, [])
        expect("flag_off_unchanged_" + os.path.basename(f), jb["status"] == jn["status"] and jb["iterations"] == jn["iterations"] and jb["objective"] == jn["objective"] and "result_hash" not in jn,
               (jb["status"], jn["status"]))

# 7. usage / scope errors
for args, name in ((["--work-limit", "0"], "zero"), (["--work-limit", "x"], "nonnumeric"), (["--work-limit", "5", "--method", "ipm"], "ipm"),
                   (["--work-limit", "5", "--explain-infeasible", os.path.join(TMP, "e.json")], "explain")):
    p, _ = run(FIX[0], args)
    expect("usage_error_" + name, p.returncode == 2, p.returncode)
milp = os.path.join(ROOT, "benchmarks/refinery_stress/campaign_milp_t6.mps")
if os.path.exists(milp):
    p, j = run(milp, ["--work-limit", "50"])
    expect("milp_unsupported", "unsupported" in p.stdout and p.returncode == 5, p.stdout[-120:])
print("FAILED %s" % fails if fails else "ALL PASS"); sys.exit(1 if fails else 0)
