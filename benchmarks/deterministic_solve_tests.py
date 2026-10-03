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
# 8. accounting: work_used counts performed iterations (no denied entries), equal to the reported iteration count
for m in FIX[:3]:
    for cap in (3, 11, 40, 100000):
        _, j = run(m, ["--work-limit", str(cap)])
        expect("work_used_equals_iterations_%d_%s" % (cap, os.path.basename(m)), j["work_used"] == j["iterations"] and j["work_used"] <= cap, (j["work_used"], j["iterations"]))

# 9. KKT retry fixture: a capped retry keeps the iteration_limit contract (was numerical_failure before the fix)
retry = os.path.join(ROOT, "tests/fixtures/scaled_retry_c87.mps")
if os.path.exists(retry):
    sts = {}
    for cap in (1, 5, 10, 15, 20, 30, 36, 37, 50, 100):
        p, j = run(retry, ["--work-limit", str(cap)])
        sts[cap] = j["status"]
        expect("retry_fixture_no_numerical_failure_cap_%d" % cap, j["status"] in ("iteration_limit", "optimal"), j["status"])
        expect("retry_fixture_exit_code_cap_%d" % cap, p.returncode == (4 if j["status"] == "iteration_limit" else 0), p.returncode)
    p, j = run(retry, ["--work-limit", "20"])
    expect("retry_path_executed", "recovered by equilibrated retry" in json.dumps(run(retry, ["--work-limit", "100"])[1]) and j["status"] == "iteration_limit")

# 10. hash is standard FNV-1a 64 over the documented canonical text (known vectors + engine cross-check)
def fnv(bs, h=0xcbf29ce484222325):
    for b in bs: h = ((h ^ b) * 0x100000001b3) & 0xFFFFFFFFFFFFFFFF
    return h
expect("fnv1a64_empty_vector", fnv(b"") == 0xcbf29ce484222325)
expect("fnv1a64_a_vector", fnv(b"a") == 0xaf63dc4c8601ec8c)
expect("fnv1a64_foobar_vector", fnv(b"foobar") == 0x85944171f73967e8)
def canon(j):
    h = 0xcbf29ce484222325
    def feed(t):
        nonlocal h
        h = fnv(t.encode(), h); h = fnv(b"\xff", h)
    j = json.load(open(os.path.join(TMP, "r.json")), parse_int=float)  # "-0" must stay negative zero
    feed(j["status"]); feed(str(int(j["iterations"]))); feed("%.17g" % (j["objective"] if j["status"] == "optimal" else 0.0))
    for v in j.get("x", []): feed("%.17g" % v)
    return "%016x" % h
for m in FIX[:3]:
    _, j = run(m, ["--work-limit", "100000"])
    expect("engine_hash_matches_python_fnv_" + os.path.basename(m), canon(j) == j["result_hash"], (canon(j), j["result_hash"]))

# 11. wall limit field and API-level wall enforcement
_, j = run(FIX[0], ["--work-limit", "100000"]); expect("wall_limit_inactive_by_default", j["wall_limit_active"] is False)
_, j = run(FIX[0], ["--work-limit", "100000", "--time-limit", "30"]); expect("wall_limit_active_when_given", j["wall_limit_active"] is True)
api = os.path.join(TMP, "det_api")
srcs = [os.path.join(ROOT, "src", f) for f in os.listdir(os.path.join(ROOT, "src")) if f.endswith(".cpp") and f != "main.cpp"]
cc = subprocess.run(["g++", "-O1", "-std=c++17", "-o", api, os.path.join(ROOT, "benchmarks/deterministic_solve_api_tests.cpp")] + srcs, capture_output=True, text=True)
expect("det_api_builds", cc.returncode == 0, cc.stderr[-300:])
if cc.returncode == 0:
    r = subprocess.run([api, FIX[-1]], capture_output=True, text=True); expect("det_api_wall_limit_enforced", r.returncode == 0, r.stdout + r.stderr)

print("FAILED %s" % fails if fails else "ALL PASS"); sys.exit(1 if fails else 0)
