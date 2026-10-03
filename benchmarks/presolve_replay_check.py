#!/usr/bin/env python3
"""Independent replay of a `taral --presolve-log` file against the ORIGINAL model.

The model is read with highspy (not the engine's reader). Each logged reduction is re-derived from the
replayer's own state; any precondition or logged number that does not match makes the replay fail.
For a log that claims infeasibility, the final replayed state must contain a violated column, empty row
or activity-range condition. Finally the kept rows/columns and the reduced bounds/constant must equal the
log. HiGHS is a reader/oracle only. Tolerance: relative 1e-9 (same constant as the engine) for the
derived conditions, 1e-9 relative for comparing logged numbers.
Usage: presolve_replay_check.py MODEL.mps LOG.json   (exit 0 pass, 1 fail)
"""
import json, math, sys
import numpy as np, highspy
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
import infeasibility_explanation_check as ic  # shared MPS duplicate-sum normalization (free format only)

REL = 1e-9
INF = math.inf

def sc(b): return 1 + (abs(b) if math.isfinite(b) else 0)
SNAP = 1e-6
def int_up(l): r = round(l); return float(r) if abs(l - r) <= SNAP else float(math.ceil(l))
def int_dn(u): r = round(u); return float(r) if abs(u - r) <= SNAP else float(math.floor(u))
def close(a, b): 
    if not (math.isfinite(a) and math.isfinite(b)): return a == b
    return abs(a - b) <= 1e-9 * (1 + abs(a) + abs(b))
def num(v): return {"inf": INF, "-inf": -INF}.get(v, v) if isinstance(v, str) else (math.nan if v is None else v)

def load(path):
    path = ic.normalize_mps(path)  # engine semantics: duplicate coefficients are summed (HiGHS keeps the first)
    h = highspy.Highs(); h.setOptionValue("output_flag", False)
    if h.readModel(path) != highspy.HighsStatus.kOk: raise RuntimeError("oracle could not read the model")
    lp = h.getLp(); n, m = lp.num_col_, lp.num_row_
    cols = []
    for j in range(n):
        s, e = lp.a_matrix_.start_[j], (lp.a_matrix_.start_[j + 1] if j + 1 < n else len(lp.a_matrix_.index_))
        cols.append([(lp.a_matrix_.index_[k], lp.a_matrix_.value_[k]) for k in range(s, e) if lp.a_matrix_.value_[k] != 0])
    isint = [bool(x == highspy.HighsVarType.kInteger) for x in lp.integrality_] if len(lp.integrality_) else [False] * n
    sense = 1 if lp.sense_ == highspy.ObjSense.kMinimize else -1
    return dict(n=n, m=m, cols=cols, cost=list(lp.col_cost_), lo=list(lp.col_lower_), up=list(lp.col_upper_),
                rlo=list(lp.row_lower_), rup=list(lp.row_upper_), isint=isint, offset=lp.offset_, maximize=sense < 0)

ROW_OPS = {"empty_row": (True, False), "singleton_row": (True, True), "redundant_row": (True, False)}
COL_OPS = {"round_int_bounds": (False, True), "fix_col": (False, True), "fix_empty_col": (False, True)}

def validate_schema(M, log):
    """Strict structure/finite/index validation before any replay. Returns a list of errors."""
    errs = []
    n, m = M["n"], M["m"]
    isnum = lambda v: isinstance(v, (int, float)) and not isinstance(v, bool)
    def okval(v):  # finite number or the strings "inf"/"-inf"; NaN/null never
        return v in ("inf", "-inf") or (isnum(v) and math.isfinite(v))
    for k in ("original", "infeasible", "ops", "kept_rows", "kept_cols", "reduced"):
        if k not in log: errs.append("missing key " + k)
    if errs: return errs
    if log["original"] != {"rows": m, "cols": n}: errs.append("original dimensions differ from the model")
    for k, op in enumerate(log["ops"]):
        t = op.get("type"); tag = "op %d %s" % (k, t)
        if t in ROW_OPS: need_row, need_col = ROW_OPS[t]
        elif t in COL_OPS: need_row, need_col = COL_OPS[t]
        else: errs.append(tag + ": unknown type"); continue
        r, c = op.get("row"), op.get("col")
        if not (isinstance(r, int) and not isinstance(r, bool) and isinstance(c, int) and not isinstance(c, bool)): errs.append(tag + ": non-integer index"); continue
        if need_row and not 0 <= r < m: errs.append(tag + ": row index %r out of range" % r)
        if not need_row and r != -1: errs.append(tag + ": unexpected row index %r" % r)
        if need_col and not 0 <= c < n: errs.append(tag + ": column index %r out of range" % c)
        if not need_col and c != -1: errs.append(tag + ": unexpected column index %r" % c)
        if not (isnum(op.get("a")) and math.isfinite(op["a"])): errs.append(tag + ": coefficient not finite")
        for key in ("v1", "v2", "v3", "v4"):
            if not okval(op.get(key)): errs.append(tag + ": %s is not a finite number or +-inf" % key)
    for key, size in (("kept_rows", m), ("kept_cols", n)):
        v = log[key]
        if not (isinstance(v, list) and all(isinstance(i, int) and 0 <= i < size for i in v) and v == sorted(set(v))): errs.append(key + ": invalid or unsorted indices")
    R = log["reduced"]
    for key in ("col_lo", "col_up", "row_lo", "row_up"):
        if not (isinstance(R.get(key), list) and all(okval(x) for x in R[key])): errs.append("reduced." + key + " invalid")
    if not okval(R.get("obj_const")): errs.append("reduced.obj_const invalid")
    if "audit" in log:
        a = log["audit"]
        for key in ("max_row_violation", "max_bound_violation", "max_int_violation", "objective", "reported_objective", "objective_diff"):
            if not (isnum(a.get(key)) and math.isfinite(a[key])): errs.append("audit." + key + " is not a finite number")
        if not isinstance(a.get("ok"), bool): errs.append("audit.ok is not a boolean")
        elif a["ok"] and not errs:
            if max(a["max_row_violation"], a["max_bound_violation"], a["max_int_violation"]) > 1e-6: errs.append("audit.ok is true but a violation exceeds 1e-6")
            if a["objective_diff"] > 1e-6 * (1 + abs(a["objective"])): errs.append("audit.ok is true but the objective differs from the reported one")
            if abs(abs(a["objective"] - a["reported_objective"]) - a["objective_diff"]) > 1e-9 * (1 + abs(a["objective"])): errs.append("audit.objective_diff is inconsistent with objective and reported_objective")
    return errs

def replay(M, log):
    pre = validate_schema(M, log)
    if pre: return pre
    return replay_ops(M, log)

def replay_ops(M, log):
    n, m = M["n"], M["m"]
    lo, up, rlo, rup = M["lo"][:], M["up"][:], M["rlo"][:], M["rup"][:]
    objc = M["offset"]
    rows = [[] for _ in range(m)]
    for j, col in enumerate(M["cols"]):
        for i, v in col: rows[i].append((j, v))
    rrem, crem = [False] * m, [False] * n
    errs = []
    def act_cols(i): return [(j, v) for j, v in rows[i] if not crem[j]]
    for k, op in enumerate(log["ops"]):
        t, i, j = op["type"], op["row"], op["col"]
        v = [num(op[x]) for x in ("v1", "v2", "v3", "v4")]
        tag = "op %d %s" % (k, t)
        if t == "round_int_bounds":
            if not M["isint"][j]: errs.append(tag + ": column not integer"); continue
            if not (close(v[0], lo[j]) and close(v[1], up[j])): errs.append(tag + ": old bounds mismatch")
            nl = int_up(lo[j]) if math.isfinite(lo[j]) else lo[j]
            nu = int_dn(up[j]) if math.isfinite(up[j]) else up[j]
            if not (close(v[2], nl) and close(v[3], nu)): errs.append(tag + ": new bounds mismatch")
            lo[j], up[j] = nl, nu
        elif t == "fix_col":
            if crem[j] or not (lo[j] == up[j] and math.isfinite(lo[j])) or not close(v[0], lo[j]): errs.append(tag + ": not fixed at logged value"); continue
            crem[j] = True
            for r, a in M["cols"][j]:
                if rrem[r]: continue
                if math.isfinite(rlo[r]): rlo[r] -= a * v[0]
                if math.isfinite(rup[r]): rup[r] -= a * v[0]
            objc += M["cost"][j] * v[0]
        elif t == "empty_row":
            if rrem[i] or act_cols(i): errs.append(tag + ": row not empty/active"); continue
            if not (close(v[0], rlo[i]) and close(v[1], rup[i])): errs.append(tag + ": bounds mismatch")
            if rlo[i] > REL * sc(rlo[i]) or rup[i] < -REL * sc(rup[i]): errs.append(tag + ": row not satisfied by 0")
            rrem[i] = True
        elif t == "singleton_row":
            ac = act_cols(i)
            if rrem[i] or len(ac) != 1 or ac[0][0] != j or not close(op["a"], ac[0][1]): errs.append(tag + ": not a singleton on logged column/coefficient"); continue
            if not (close(v[0], rlo[i]) and close(v[1], rup[i])): errs.append(tag + ": row bounds mismatch")
            a = ac[0][1]
            l, u = (rlo[i] / a, rup[i] / a) if a > 0 else (rup[i] / a, rlo[i] / a)
            if M["isint"][j]:
                if math.isfinite(l): l = int_up(l)
                if math.isfinite(u): u = int_dn(u)
            nl, nu = max(lo[j], l), min(up[j], u)
            if (nl > nu) if M["isint"][j] else (nl > nu + REL * sc(nu)): errs.append(tag + ": conflict logged as reduction")
            if nl > nu: nu = nl
            if not (close(v[2], nl) and close(v[3], nu)): errs.append(tag + ": derived bounds mismatch")
            lo[j], up[j] = nl, nu; rrem[i] = True
        elif t == "redundant_row":
            if rrem[i]: errs.append(tag + ": row already removed"); continue
            mn = sum((a * lo[c] if a > 0 else a * up[c]) for c, a in act_cols(i))
            mx = sum((a * up[c] if a > 0 else a * lo[c]) for c, a in act_cols(i))
            if not (close(v[0], rlo[i]) and close(v[1], rup[i]) and close(v[2], mn) and close(v[3], mx)): errs.append(tag + ": activity/bounds mismatch")
            ok = (math.isfinite(mn) and math.isfinite(mx) and mn >= rlo[i] and mx <= rup[i]) or \
                 (rlo[i] == -INF and math.isfinite(mx) and mx <= rup[i]) or (rup[i] == INF and math.isfinite(mn) and mn >= rlo[i])
            if not ok: errs.append(tag + ": row is not implied by column bounds")
            rrem[i] = True
        elif t == "fix_empty_col":
            if crem[j] or any(not rrem[r] for r, a in M["cols"][j]): errs.append(tag + ": column still in an active row"); continue
            c = -M["cost"][j] if M["maximize"] else M["cost"][j]
            want = lo[j] if c > 0 else (up[j] if c < 0 else (lo[j] if math.isfinite(lo[j]) else (up[j] if math.isfinite(up[j]) else 0.0)))
            if not math.isfinite(want) or not close(v[0], want): errs.append(tag + ": value is not the cost-optimal bound")
            lo[j] = up[j] = v[0]
        else:
            errs.append(tag + ": unknown op")
    if log.get("timed_out"):
        return errs  # incomplete presolve: operations validated, no reduced model is claimed
    if log["infeasible"]:
        found = any((lo[c] > up[c]) if M["isint"][c] else (lo[c] > up[c] + REL * sc(up[c])) for c in range(n) if not crem[c])
        for i in range(m):
            if rrem[i] or found: continue
            ac = act_cols(i)
            if not ac:
                found = rlo[i] > REL * sc(rlo[i]) or rup[i] < -REL * sc(rup[i])
            else:
                mn = sum((a * lo[c] if a > 0 else a * up[c]) for c, a in ac); mx = sum((a * up[c] if a > 0 else a * lo[c]) for c, a in ac)
                found = (math.isfinite(mn) and mn > rup[i] + REL * sc(rup[i])) or (math.isfinite(mx) and mx < rlo[i] - REL * sc(rlo[i]))
                if len(ac) == 1:  # singleton conflict against bounds
                    a = ac[0][1]; l, u = (rlo[i] / a, rup[i] / a) if a > 0 else (rup[i] / a, rlo[i] / a)
                    c = ac[0][0]
                    if M["isint"][c]:
                        if math.isfinite(l): l = int_up(l)
                        if math.isfinite(u): u = int_dn(u)
                        found = found or max(lo[c], l) > min(up[c], u)
                    else:
                        found = found or max(lo[c], l) > min(up[c], u) + REL * sc(min(up[c], u))
        if not found: errs.append("log claims infeasible but the replayed state shows no violated condition")
        return errs
    kr = [i for i in range(m) if not rrem[i]]; kc = [c for c in range(n) if not crem[c]]
    if kr != log["kept_rows"] or kc != log["kept_cols"]: errs.append("kept rows/cols differ from the log")
    R = log["reduced"]
    if not close(R["obj_const"], objc): errs.append("reduced objective constant differs")
    for name, mine, idx in (("col_lo", lo, kc), ("col_up", up, kc), ("row_lo", rlo, kr), ("row_up", rup, kr)):
        got = [num(x) for x in R[name]]
        if len(got) != len(idx) or any(not close(got[q], mine[idx[q]]) for q in range(len(idx))): errs.append("reduced %s differs" % name)
    return errs

if __name__ == "__main__":
    def bad_const(c): raise ValueError("non-standard JSON constant " + c)
    try:
        M = load(sys.argv[1]); log = json.load(open(sys.argv[2]), parse_constant=bad_const)
        errs = replay(M, log)
        if not errs and len(sys.argv) > 3 and "audit" in log:  # optional: recompute the audited objective from the written .sol
            val = dict(l.split() for l in open(sys.argv[3]))
            h = highspy.Highs(); h.setOptionValue("output_flag", False); h.readModel(ic.normalize_mps(sys.argv[1])); names = list(h.getLp().col_names_)
            x = [float(val[nm]) for nm in names]
            obj = M["offset"] + sum(c * v for c, v in zip(M["cost"], x))
            if abs(obj - log["audit"]["objective"]) > 1e-6 * (1 + abs(obj)): errs.append("audit objective %.9g != objective recomputed from the sol file %.9g" % (log["audit"]["objective"], obj))
    except Exception as ex:
        print("replay FAIL\n  - %s: %s" % (type(ex).__name__, ex)); sys.exit(1)
    print("replay %s (%d ops)" % ("PASS" if not errs else "FAIL", len(log["ops"])))
    for e in errs[:10]: print("  -", e)
    sys.exit(1 if errs else 0)
