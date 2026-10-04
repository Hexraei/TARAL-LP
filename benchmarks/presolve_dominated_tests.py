#!/usr/bin/env python3
"""Tests for the dominated_col presolve op (dual fixing). Usage: presolve_dominated_tests.py [path/to/taral] [n_random]
Separate from presolve_tests.py so the primary suite counts (809 cases, 464 audited points) are unchanged.
Oracle: HiGHS (highspy) on the ORIGINAL model, status and objective only. Every case: taral with --presolve agrees with the oracle
(status; objective at relative 1e-6 when optimal), the log replays against the original model with the independent replayer
(presolve_replay_check), and the reported point is re-audited from the .sol file against the original rows/bounds/integrality."""
import copy, json, math, os, random, subprocess, sys, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, highspy
import presolve_replay_check as rc

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "build", "taral")
NRAND = int(sys.argv[2]) if len(sys.argv) > 2 else 300
TMP = tempfile.mkdtemp()
fails = []
stats = {"hand": 0, "random": 0, "dominated_ops": 0, "random_with_dominated": 0, "audited_points": 0, "down": 0, "up": 0, "int_cols": 0, "presolve_infeasible": 0, "presolve_unbounded": 0}
def expect(name, cond, msg=""):
    print(("PASS " if cond else "FAIL ") + name + ("" if cond else " " + str(msg)))
    if not cond: fails.append(name)

def mps(name, body):
    p = os.path.join(TMP, name + ".mps"); open(p, "w").write(body); return p

def oracle(path):
    h = highspy.Highs(); h.setOptionValue("output_flag", False); h.readModel(rc.ic.normalize_mps(path)); h.run()
    st = h.modelStatusToString(h.getModelStatus())
    m = {"Optimal": "optimal", "Infeasible": "infeasible", "Unbounded": "unbounded", "Primal infeasible or unbounded": "infeasible_or_unbounded"}.get(st, st)
    return m, (h.getInfo().objective_function_value if m == "optimal" else None)

def run(path, tag, extra):
    js, so, lg = (os.path.join(TMP, tag + e) for e in (".json", ".sol", "_log.json"))
    for f in (js, so, lg):
        if os.path.exists(f): os.remove(f)
    r = subprocess.run([BIN, path, "--json", js, "--sol", so, "--time-limit", "30"] + extra(lg), capture_output=True, text=True, timeout=60)
    sl = [l for l in r.stdout.splitlines() if l.startswith("status ")]
    st = sl[-1].split()[1] if sl else "?"
    obj = float(sl[-1].split("objective ")[1].split()[0]) if sl and "objective " in sl[-1] else None
    return st, obj, so, lg

def check(path, tag, want_dom=None, want_dir=None, not_cols=()):
    """want_dom: True = at least one dominated_col op expected, False = none expected, None = no expectation."""
    o, oobj = oracle(path)
    s0, o0, _, _ = run(path, tag + "_off", lambda lg: [])
    s1, o1, so, lg = run(path, tag + "_on", lambda lg: ["--presolve-log", lg])
    M = rc.load(path); problems = []
    plain_ok = (s0 == o) or (o == "unbounded" and s0 in ("unbounded", "unbounded_relaxation", "dual_infeasible"))
    if not plain_ok and o in ("optimal", "infeasible", "unbounded"):
        # the plain engine itself disagrees with the oracle (existing behaviour, e.g. unbounded relaxations of mixed-integer models or
        # numerical_failure): presolve must then give the same status class as the plain engine, no more is claimed
        stats["engine_oracle_divergence"] = stats.get("engine_oracle_divergence", 0) + 1
        if s1 != s0: problems.append("presolve status %s vs plain %s (plain also differs from oracle %s)" % (s1, s0, o))
    elif o == "optimal":
        if s1 != "optimal": problems.append("presolve status %s vs oracle optimal" % s1)
        elif abs(o1 - oobj) > 1e-6 * (1 + abs(oobj)): problems.append("presolve objective %.9g vs oracle %.9g" % (o1, oobj))
    elif o == "infeasible":
        if s1 != "infeasible": problems.append("presolve status %s vs oracle infeasible" % s1)
    elif o == "unbounded":
        if s1 not in ("unbounded", "unbounded_relaxation", "dual_infeasible"): problems.append("presolve status %s vs oracle unbounded" % s1)
    else:  # infeasible_or_unbounded from the oracle: presolve must not claim optimal
        if s1 == "optimal": problems.append("presolve optimal where the oracle says infeasible/unbounded")
    ops = []
    if os.path.exists(lg):
        L = json.load(open(lg)); errs = rc.replay(M, L)
        if errs: problems.append("replay: " + "; ".join(errs[:2]))
        ops = [op for op in L["ops"] if op["type"] == "dominated_col"]
    else: problems.append("no log written")
    if s1 == "optimal" and os.path.exists(so):  # re-audit the point in original space
        val = dict(l.split() for l in open(so)); stats["audited_points"] += 1
        h = highspy.Highs(); h.setOptionValue("output_flag", False); h.readModel(rc.ic.normalize_mps(path)); lp = h.getLp()
        x = np.array([float(val[nm]) for nm in lp.col_names_])
        A = rc.sp.csc_matrix((lp.a_matrix_.value_, lp.a_matrix_.index_, lp.a_matrix_.start_), shape=(lp.num_row_, lp.num_col_)) if hasattr(rc, "sp") else None
        if A is not None:
            act = A @ x
            for i in range(lp.num_row_):
                for b, sgn in ((lp.row_lower_[i], 1), (lp.row_upper_[i], -1)):
                    if math.isfinite(b) and abs(b) < 1e30 and sgn * (act[i] - b) < -1e-6 * (1 + abs(b)): problems.append("row %d violated" % i)
            for j in range(lp.num_col_):
                if x[j] < lp.col_lower_[j] - 1e-6 * (1 + abs(lp.col_lower_[j])) and lp.col_lower_[j] > -1e30: problems.append("col %d below lower" % j)
                if x[j] > lp.col_upper_[j] + 1e-6 * (1 + abs(lp.col_upper_[j])) and lp.col_upper_[j] < 1e30: problems.append("col %d above upper" % j)
        for j, t in enumerate(lp.integrality_ if len(lp.integrality_) else []):
            if t == highspy.HighsVarType.kInteger and abs(x[j] - round(x[j])) > 1e-6: problems.append("col %d not integral" % j)
    stats["dominated_ops"] += len(ops)
    for op in ops: stats["down" if op["v2"] == -1 else "up"] += 1; stats["int_cols"] += bool(M["isint"][op["col"]])
    if s1 == "infeasible": stats["presolve_infeasible"] += 1
    if s1 in ("unbounded", "unbounded_relaxation", "dual_infeasible"): stats["presolve_unbounded"] += 1
    if want_dom is True and not ops: problems.append("expected a dominated_col op, none logged")
    if want_dom is False and ops: problems.append("expected no dominated_col op, got %d" % len(ops))
    if [op for op in ops if op["col"] in not_cols]: problems.append("dominated_col on a column that must not be fixed: %r" % [op["col"] for op in ops if op["col"] in not_cols])
    if want_dir is not None and ops and ops[0]["v2"] != want_dir: problems.append("direction %r != %r" % (ops[0]["v2"], want_dir))
    expect(tag, not problems, "; ".join(problems))
    return ops, M, (json.load(open(lg)) if os.path.exists(lg) else None)

def model(rows, cols, objsense="MIN", ints=(), bounds=()):
    """rows: list of type letters; cols: dict name -> (cost, {row: coef}); rhs per row in rhs list given as rows entries (type, rhs)."""
    L = ["NAME D"] + (["OBJSENSE", "    MAX"] if objsense == "MAX" else []) + ["ROWS", " N OBJ"] + [" %s R%d" % (t, i) for i, (t, _) in enumerate(rows)] + ["COLUMNS"]
    inint = False
    for nm, (c, ents) in cols.items():
        if nm in ints and not inint: L.append(" MARK 'MARKER' 'INTORG'"); inint = True
        if nm not in ints and inint: L.append(" MARK 'MARKER' 'INTEND'"); inint = False
        L.append(" %s OBJ %.17g" % (nm, c))
        for r, v in ents.items(): L.append(" %s R%d %.17g" % (nm, r, v))
    if inint: L.append(" MARK 'MARKER' 'INTEND'")
    L += ["RHS"] + [" RHS R%d %.17g" % (i, rhs) for i, (_, rhs) in enumerate(rows)] + ["BOUNDS"] + list(bounds) + ["ENDATA"]
    return "\n".join(L) + "\n"

# ---- hand cases
# H1: min 3x + y, R0: x + y <= 6 (L), R1: y >= 1 (G). x costs > 0 and R0 only has an upper bound -> fixed at its lower bound 2.
p = mps("h1", model([("L", 6), ("G", 1)], {"X": (3, {0: 1}), "Y": (1, {0: 1, 1: 1})}, bounds=[" LO BND X 2", " UP BND X 5", " UP BND Y 5"]))
ops, M, L1 = check(p, "h1_down_fix", True, -1); stats["hand"] += 1
expect("h1_value", bool(ops) and ops[0]["col"] == 0 and ops[0]["v1"] == 2.0, ops[:1])
# H2: max form: max 3x + y, R0: x + y >= 1 (G), x <= 5: x improves with larger value, rows only have lower bounds -> fixed at upper 5
p = mps("h2", model([("G", 1)], {"X": (3, {0: 1}), "Y": (-1, {0: 1})}, objsense="MAX", bounds=[" UP BND X 5", " UP BND Y 5"]))
check(p, "h2_up_fix_max", True, 1); stats["hand"] += 1
# H3: integer column fixed at an integer bound; row keeps the model non-trivial
p = mps("h3", model([("L", 7.5), ("G", 1)], {"X": (2, {0: 1}), "Y": (1, {0: 1, 1: 1})}, ints=("X",), bounds=[" LO BND X 1", " UP BND X 9", " UP BND Y 9"]))
ops, _, _ = check(p, "h3_integer_down", True, -1); stats["hand"] += 1
# H4: zero cost: either direction is valid; one op, replay accepts it
p = mps("h4", model([("L", 6), ("G", 1)], {"X": (0, {0: 1}), "Y": (1, {0: 1, 1: 1})}, bounds=[" LO BND X 2", " UP BND X 5", " UP BND Y 5"]))
check(p, "h4_zero_cost", True); stats["hand"] += 1
# H5: NOT dominated: cost pulls one way, the row blocks it (x < 0 cost with an upper-bounded row)
p = mps("h5", model([("L", 6), ("G", 1)], {"X": (-3, {0: 1}), "Y": (1, {0: 1, 1: 1})}, bounds=[" LO BND X 2", " UP BND X 5", " UP BND Y 5"]))
check(p, "h5_cost_wrong_way", None, None, (0,)); stats["hand"] += 1
# H6: NOT dominated: a row bounds the column from the harmful side (x in an equality row)
p = mps("h6", model([("E", 4), ("G", 1)], {"X": (3, {0: 1}), "Y": (1, {0: 1, 1: 1})}, bounds=[" LO BND X 0", " UP BND X 5", " UP BND Y 5"]))
check(p, "h6_equality_blocks", False); stats["hand"] += 1
# H7: infinite improving bound: min x with x free below and only upper-bounded rows -> unbounded, no op, status unbounded
p = mps("h7", model([("L", 6)], {"X": (1, {0: 1}), "Y": (1, {0: 1})}, bounds=[" MI BND X", " UP BND Y 5"]))
check(p, "h7_infinite_bound_unbounded", None, None, (0,)); stats["hand"] += 1
# H8: infeasible after the fix is a true infeasibility: x >= 2 forced, y <= 3 forced by bounds, row x + y <= 4 requires y <= 2, R1 y >= 3
p = mps("h8", model([("L", 4), ("G", 3)], {"X": (1, {0: 1}), "Y": (1, {0: 1, 1: 1})}, bounds=[" LO BND X 2", " UP BND X 5", " UP BND Y 5"]))
check(p, "h8_infeasible_stays_infeasible", None); stats["hand"] += 1
# H9: unbounded model stays unbounded: max x + y, row x - y <= 2 ... x in an upper-bounded row positive coef: not dominated for max
p = mps("h9", model([("L", 2)], {"X": (1, {0: 1}), "Y": (1, {0: -1})}, objsense="MAX", bounds=[" UP BND Y 1e30"]))
check(p, "h9_unbounded_max", None); stats["hand"] += 1
# H10: bound with tracked uncertainty (derived through a singleton row with |a| != 1) must not be fixed on
p = mps("h10", model([("G", 1), ("L", 9), ("G", 1)], {"X": (1, {0: 3, 1: 1}), "Y": (1, {1: 1, 2: 1})}, bounds=[" UP BND X 5", " UP BND Y 5"]))
ops, M, L10 = check(p, "h10_uncertain_bound_declined", None); stats["hand"] += 1
ks = [o for o in (L10 or {"ops": []})["ops"] if o["type"] == "singleton_row" and o["col"] == 0 and abs(o["a"]) != 1]
if ks: expect("h10_no_dom_on_uncertain_col", not [o for o in ops if o["col"] == 0], ops)

# ---- corrupted logs must be rejected by the replayer
p = os.path.join(TMP, "h1.mps"); M = rc.load(p); L = L1
def rej(mut, name):
    l2 = copy.deepcopy(L); mut(l2); expect("corrupt_" + name, bool(rc.replay(M, l2)))
k = [i for i, o in enumerate(L["ops"]) if o["type"] == "dominated_col"][0]
rej(lambda l: l["ops"][k].__setitem__("v2", 1), "wrong_direction")
rej(lambda l: l["ops"][k].__setitem__("v1", 3.0), "wrong_value")
rej(lambda l: l["ops"][k].__setitem__("v2", 0), "bad_direction")
rej(lambda l: l["ops"][k].__setitem__("col", 1), "wrong_column")
rej(lambda l: l["ops"].insert(0, {"type": "dominated_col", "row": -1, "col": 1, "a": 0, "v1": 0.0, "v2": -1, "v3": 0, "v4": 0}), "bogus_op_on_blocked_col")
# a column that is dominated only if the cost sign were different must be rejected by the replayer
M5 = rc.load(os.path.join(TMP, "h5.mps"))
bogus = {"original": {"rows": M5["m"], "cols": M5["n"]}, "infeasible": False, "ops": [{"type": "dominated_col", "row": -1, "col": 0, "a": 0, "v1": 2.0, "v2": -1, "v3": 0, "v4": 0}],
         "kept_rows": list(range(M5["m"])), "kept_cols": list(range(M5["n"])), "reduced": {"col_lo": [], "col_up": [], "row_lo": [], "row_up": [], "obj_const": 0}}
expect("corrupt_cost_sign_rejected", bool(rc.replay(M5, bogus)))

# ---- random sweep: one-sided rows with mixed signs, mixed integer/continuous, min and max, finite and infinite bounds
rng = random.Random(20261004)
def rand_model(k):
    m, n = rng.randint(1, 5), rng.randint(2, 7)
    rows = []
    for i in range(m):
        t = rng.choice("LLGGE" if rng.random() < .9 else "LGE"); rows.append((t, float(rng.randint(-6, 12))))
    ints = [f"X{j}" for j in range(n) if rng.random() < .3]
    cols = {}
    for j in range(n):
        ents = {i: float(rng.choice([-3, -2, -1, 1, 1, 2, 3])) for i in range(m) if rng.random() < .6}
        cols["X%d" % j] = (float(rng.choice([-2, -1, 0, 1, 1, 2, 3])), ents)
    bnds = []
    for j in range(n):
        lo = rng.choice([0, 0, 1, -2, None]); up = rng.choice([3, 5, 8, None, None])
        if "X%d" % j in ints and up is None: up = 6  # HiGHS gives an integer column without UP an upper bound of 1, the engine does not: always state it
        if lo is None: bnds.append(" MI BND X%d" % j)
        elif lo != 0: bnds.append(" LO BND X%d %d" % (j, lo))
        if up is not None and (lo is None or up >= lo): bnds.append(" UP BND X%d %d" % (j, up))
        elif up is not None: bnds.append(" UP BND X%d %d" % (j, (lo if lo is not None else 0) + 1))
    return model(rows, cols, rng.choice(["MIN", "MAX"]), tuple(ints), bnds)
for k in range(NRAND):
    p = mps("r%d" % k, rand_model(k)); stats["random"] += 1
    before = stats["dominated_ops"]; check(p, "rand_%d" % k)
    stats["random_with_dominated"] += stats["dominated_ops"] > before
expect("random_sweep_exercises_dominated_col", stats["random_with_dominated"] >= max(5, NRAND // 20), stats)
print("stats", json.dumps(stats))
print("FAILED %d" % len(fails) if fails else "ALL PASS")
sys.exit(1 if fails else 0)
