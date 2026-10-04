#!/usr/bin/env python3
"""Independent check of `taral --explain-infeasible` output. HiGHS is a test oracle only.

Reads the MPS with highspy (not the engine's reader) and checks, from the JSON alone:
  1. status 'irreducible': the explained rows plus all column bounds are infeasible, and
     dropping any one explained row makes them feasible (oracle LP solves), and each stored
     removal witness satisfies every other explained row and all column bounds (rel tol 1e-7).
  2. the relaxation point satisfies all column bounds, and each row within its stated relaxation;
     its weighted L1 objective matches the stated one and equals the oracle optimum (rel 1e-6).
Scope is the LP relaxation (integrality ignored). 'Irreducible' = row-irreducible under retained
column bounds; it is not claimed globally minimum.
Usage: infeasibility_explanation_check.py MODEL.mps EXPLANATION.json   (exit 0 pass, 1 fail)
"""
import json, math, os, sys, tempfile
import numpy as np, highspy, scipy.sparse as sp
from scipy.optimize import linprog

TOL = 1e-7

def normalize_mps(path):
    """Sum duplicate (column,row) coefficients in COLUMNS (the engine's semantics; highspy keeps the first).
    Also uses only the first BOUNDS set and drops an earlier entry only when a later entry of the same column sets every side it sets (the engine overwrites; HiGHS keeps the first and warns); entries are never rewritten, and a kept overlap raises. Free-format token lines only; MARKER lines and other sections are copied unchanged."""
    out, sec, seen, cur, bent, bset = [], None, {}, None, [], None
    for raw in open(path).read().splitlines():
        if not raw.strip() or raw.startswith("*"): out.append(raw); continue
        if not raw[0].isspace():
            if sec == "COLUMNS": out.extend(flush(seen))
            sec = raw.split()[0]; seen = {}; out.append(raw); continue
        if sec == "BOUNDS":  # engine semantics: only the first bound set is used; entries are kept IN ORDER (order-dependent rules such as a negative UP stay with HiGHS/engine)
            t = raw.split(); ty = t[0].upper()
            noval = ty in ("MI", "PL", "FR")
            if ty in ("LO", "UP", "FX", "MI", "PL", "FR") and len(t) in ((2, 3) if noval else (3, 4)):
                has_set = len(t) == (3 if noval else 4)
                if has_set:
                    if bset is None: bset = t[1]
                    if t[1] != bset: continue
                bent.append((len(out), t[2] if has_set else t[1], ty)); out.append(raw)
                continue
            out.append(raw); continue
        if sec != "COLUMNS": out.append(raw); continue
        t = raw.split()
        if "'MARKER'" in t or "MARKER" in t: out.extend(flush(seen)); seen = {}; out.append(raw); continue
        if len(t) not in (3, 5): raise ValueError("unsupported COLUMNS line (need free format): " + raw)
        for k in range(1, len(t), 2):
            seen.setdefault(t[0], {}); seen[t[0]][t[k]] = seen[t[0]].get(t[k], 0.0) + float(t[k + 1].replace("D", "E").replace("d", "e"))
    # Drop an entry only when every side it sets (lower/upper) is set again by a LATER entry of the same column; entries are
    # never rewritten, so order-dependent rules (negative UP with an implicit lower, ...) are left to the readers.
    sides = {"LO": {"lo"}, "MI": {"lo"}, "UP": {"up"}, "PL": {"up"}, "FX": {"lo", "up"}, "FR": {"lo", "up"}}
    later, drop, kept_sides = {}, set(), {}
    for pos, col, ty in reversed(bent):
        cov = later.setdefault(col, set())
        if sides[ty] <= cov: drop.add(pos); continue
        cov |= sides[ty]; kept_sides.setdefault(col, []).append(ty)
    for col, tys in kept_sides.items():
        flat = [x for ty in tys for x in sides[ty]]
        if len(flat) != len(set(flat)): raise ValueError("unsupported BOUNDS combination for column %s (%s): kept entries overlap" % (col, tys))
    # Engine rule for a negative UP bound: if the lower bound is still 0 at that point (default or an explicit 0), the lower becomes -inf
    # (the MPS convention); highspy keeps lower 0 and reports infeasible, so the equivalent MI entry is inserted before such an UP.
    lower, pre, lastlow = {}, {}, {}
    # An earlier negative UP can be dropped because a later UP sets the same (upper) side, but its engine effect on the LOWER bound
    # (lower -> -inf when the lower is still 0) survives that overwrite: the engine's UP -1 then UP 7 gives [-inf, 7], not [0, 7].
    # So the pass runs over ALL entries in file order; a dropped negative UP that triggers the rule is replaced by the equivalent MI.
    later_low = {}  # does a later entry of the column set the lower bound? then the engine's lower is overwritten anyway
    for pos, col, ty in reversed(bent):
        later_low[pos] = later_low.get(("c", col), False)
        if ty in ("LO", "FX", "MI", "FR"): later_low[("c", col)] = True
    for pos, col, ty in bent:
        was_dropped = pos in drop
        t = out[pos].split()
        val = float(t[-1].replace("D", "E").replace("d", "e")) if ty in ("LO", "UP", "FX") else None
        lo = lower.get(col, 0.0)
        if was_dropped:
            # a dropped entry never reaches the readers; only a triggering negative UP leaves a lower-side effect to preserve
            if not (ty == "UP" and val < 0 and lo == 0.0 and not later_low[pos]): continue
        if ty in ("LO", "FX"): lower[col] = val; lastlow[col] = (pos, ty)
        elif ty in ("MI", "FR"): lower[col] = -math.inf; lastlow[col] = (pos, ty)
        elif ty == "UP" and val < 0 and lo == 0.0 and not later_low[pos]:
            if col in lastlow:
                if lastlow[col][1] != "LO": raise ValueError("unsupported BOUNDS combination for column %s (FX 0 before a negative UP)" % col)
                drop.add(lastlow[col][0])  # explicit LO 0 is overwritten by the engine's -inf; a duplicate lower entry would be ignored by highspy
            pre[pos] = " MI %s %s" % (bset or "BND", col); lower[col] = -math.inf
    out2 = []
    for k, l in enumerate(out):
        if k in pre: out2.append(pre[k])  # before the drop test: a dropped negative UP still leaves its MI behind
        if k in drop: continue
        out2.append(l)
    out = out2
    fd, tmp = tempfile.mkstemp(suffix=".mps"); os.close(fd)
    open(tmp, "w").write("\n".join(out) + "\n")
    return tmp

def flush(seen):
    L = []
    for c, d in seen.items():
        for r, v in d.items(): L.append(" %s %s %.17g" % (c, r, v))
    return L

def load(path):
    path = normalize_mps(path)
    h = highspy.Highs(); h.setOptionValue("output_flag", False)
    if h.readModel(path) != highspy.HighsStatus.kOk: raise RuntimeError("oracle could not read the model")
    lp = h.getLp()
    n, m = lp.num_col_, lp.num_row_
    A = sp.csc_matrix((lp.a_matrix_.value_, lp.a_matrix_.index_, lp.a_matrix_.start_), shape=(m, n)).tocsr()
    return A, np.array(lp.row_lower_), np.array(lp.row_upper_), np.array(lp.col_lower_), np.array(lp.col_upper_), list(lp.row_names_)

def feasible(A, rl, ru, cl, cu, rows):
    sub = A[rows]
    Aub, bub = [], []
    Aeq, beq = [], []
    for k, i in enumerate(rows):
        lo, up = rl[i], ru[i]
        if lo == up: Aeq.append(k); beq.append(lo); continue
        if np.isfinite(up): Aub.append((k, 1.0)); bub.append(up)
        if np.isfinite(lo): Aub.append((k, -1.0)); bub.append(-lo)
    n = A.shape[1]
    Au = sp.vstack([s * sub[k] for (k, s) in Aub]) if Aub else None
    Ae = sub[Aeq] if Aeq else None
    r = linprog(np.zeros(n), A_ub=Au, b_ub=bub or None, A_eq=Ae, b_eq=beq or None,
                bounds=list(zip(np.where(np.isfinite(cl), cl, None), np.where(np.isfinite(cu), cu, None))), method="highs")
    if r.status not in (0, 2): raise RuntimeError("oracle LP ended with status %d (neither optimal nor infeasible)" % r.status)
    return r.status == 0, r.status

def viol(A, rl, ru, cl, cu, x, rows):
    ax = A[rows] @ x
    v = 0.0
    for k, i in enumerate(rows):
        s = 1 + max([abs(b) for b in (rl[i], ru[i]) if np.isfinite(b)] or [0])
        v = max(v, (rl[i] - ax[k]) / s, (ax[k] - ru[i]) / s)
    for j in range(len(x)):
        s = 1 + max([abs(b) for b in (cl[j], cu[j]) if np.isfinite(b)] or [0])
        v = max(v, (cl[j] - x[j]) / s, (x[j] - cu[j]) / s)
    return v

def check(mps, jpath):
    A, rl, ru, cl, cu, names = load(mps)
    def bad_const(c): raise ValueError("non-standard JSON constant " + c)
    e = json.load(open(jpath), parse_constant=bad_const)
    errs = []
    st = e["status"]
    m, n = A.shape
    fin = lambda v: isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v)
    def vec(v, size, nonneg=False):
        return isinstance(v, list) and len(v) == size and all(fin(a) and (a >= 0 or not nonneg) for a in v)
    for r in e["rows"]:
        if not (isinstance(r.get("index"), int) and 0 <= r["index"] < m): return st, 0, ["invalid row index %r" % r.get("index")]
    rows = [r["index"] for r in e["rows"]]
    if len(set(rows)) != len(rows): return st, 0, ["duplicate explained rows"]
    if st in ("irreducible", "reduced_unproven", "bounds_only"):
        # independent Farkas verification from the ORIGINAL model; the verified flag is never trusted
        yl, yu, zl, zu = (e.get("farkas_row_lower"), e.get("farkas_row_upper"), e.get("farkas_col_lower"), e.get("farkas_col_upper"))
        if not (vec(yl, m, True) and vec(yu, m, True) and vec(zl, n, True) and vec(zu, n, True)):
            return st, len(rows), ["Farkas multipliers: wrong dimension, non-finite or negative entries"]
        yl, yu, zl, zu = (np.array(v, float) for v in (yl, yu, zl, zu))
        l1 = yl.sum() + yu.sum() + zl.sum() + zu.sum()
        if not l1 > 0: errs.append("Farkas multipliers are all zero")
        else:
            for lab, mult, bound in (("row_lower", yl, rl), ("row_upper", yu, ru), ("col_lower", zl, cl), ("col_upper", zu, cu)):
                if np.any((mult > 0) & ~np.isfinite(bound)): errs.append("multiplier on an infinite bound (%s)" % lab)
            if not errs:
                stat = A.T @ (yu - yl) + zu - zl                      # must vanish: the combination has no x term
                resid = float(np.max(np.abs(stat))) / l1
                gap = (float(np.dot(yu, np.where(yu > 0, ru, 0))) - float(np.dot(yl, np.where(yl > 0, rl, 0)))
                       + float(np.dot(zu, np.where(zu > 0, cu, 0))) - float(np.dot(zl, np.where(zl > 0, cl, 0))))
                margin = -gap / l1                                     # >0: 0 = combination <= gap < 0, contradiction
                if resid > 1e-9: errs.append("Farkas stationarity residual %.3g > 1e-9" % resid)
                if not margin > 1e-8: errs.append("Farkas contradiction margin %.3g not > 1e-8" % margin)
                support = set(np.nonzero(yl + yu)[0].tolist())
                if st != "bounds_only" and not support <= set(rows): errs.append("multipliers use rows outside the explained set")
        if errs: return st, len(rows), errs
    for r in e["rows"]:
        if names and names[r["index"]] != r["name"]: errs.append("row name mismatch %s" % r["name"])
    if st == "irreducible":
        ok, code = feasible(A, rl, ru, cl, cu, rows)
        if ok: errs.append("explained rows are feasible in the oracle")
        for pos, i in enumerate(rows):
            rest = [r for r in rows if r != i]
            ok, _ = feasible(A, rl, ru, cl, cu, rest)
            if not ok: errs.append("row %s is removable (rest still infeasible)" % e["rows"][pos]["name"])
            w = e["rows"][pos]["removal_witness"]
            if not vec(w, n): errs.append("witness for %s is missing, wrong length or non-finite" % e["rows"][pos]["name"]); continue
            v = viol(A, rl, ru, cl, cu, np.array(w), rest)
            if v > TOL: errs.append("witness for %s violates by %.3g" % (e["rows"][pos]["name"], v))
    elif st == "bounds_only":
        if not e["inconsistent_bound_columns"]: errs.append("bounds_only without inconsistent columns")
        for j in e["inconsistent_bound_columns"]:
            if not cl[j] > cu[j]: errs.append("column %d bounds are consistent" % j)
    elif st == "relaxation_feasible":
        ok, _ = feasible(A, rl, ru, cl, cu, list(range(A.shape[0])))
        if not ok: errs.append("model claimed feasible but oracle infeasible")
    # relaxation (only meaningful when infeasible)
    rx = e["relaxation"]
    if rx["status"] == "optimal" and st in ("irreducible", "reduced_unproven"):
        if not (vec(rx["x"], n) and all(isinstance(r.get("index"), int) and 0 <= r["index"] < m and fin(r["lower_relaxed_by"]) and fin(r["upper_relaxed_by"])
                and r["lower_relaxed_by"] >= 0 and r["upper_relaxed_by"] >= 0 for r in rx["rows"]) and fin(rx["objective"])):
            errs.append("relaxation: invalid x/slack/objective values"); return st, len(rows), errs
        x = np.array(rx["x"])
        sl, su = np.zeros(m), np.zeros(m); w = np.zeros(m)
        for r in rx["rows"]: sl[r["index"]] = r["lower_relaxed_by"]; su[r["index"]] = r["upper_relaxed_by"]; w[r["index"]] = r["weight"]
        ax = A @ x
        for i in range(m):
            if ax[i] < rl[i] - sl[i] - TOL * (1 + abs(rl[i]) if np.isfinite(rl[i]) else 1) or ax[i] > ru[i] + su[i] + TOL * (1 + abs(ru[i]) if np.isfinite(ru[i]) else 1):
                errs.append("relaxation point outside relaxed row %d" % i); break
        for j in range(n):
            if x[j] < cl[j] - TOL * (1 + abs(cl[j])) or x[j] > cu[j] + TOL * (1 + abs(cu[j])): errs.append("relaxation point violates column bound %d" % j); break
        wt = np.array([1.0 / (1 + max([abs(b) for b in (rl[i], ru[i]) if np.isfinite(b)] or [0])) for i in range(m)])
        # recompute the minimal relaxation of the point itself, independent of the stated slacks
        pl = np.maximum(0, rl - ax); pu = np.maximum(0, ax - ru)
        pobj = float(wt @ (pl + pu))
        if abs(pobj - rx["objective"]) > 1e-6 * (1 + abs(pobj)): errs.append("objective %.9g != recomputed %.9g" % (rx["objective"], pobj))
        # oracle optimum: min sum w (s_lo + s_up)
        c = np.concatenate([np.zeros(n), wt, wt])
        I = sp.identity(m, format="csr")
        rows_ub, b_ub = [], []
        fl = np.isfinite(rl); fu = np.isfinite(ru)
        Aub = sp.vstack([sp.hstack([-A[fl], -I[fl], sp.csr_matrix((fl.sum(), m))]),
                         sp.hstack([A[fu], sp.csr_matrix((fu.sum(), m)), -I[fu]])])
        b = np.concatenate([-rl[fl], ru[fu]])
        bnds = [(cl[j] if np.isfinite(cl[j]) else None, cu[j] if np.isfinite(cu[j]) else None) for j in range(n)] + [(0, None)] * (2 * m)
        r = linprog(c, A_ub=Aub, b_ub=b, bounds=bnds, method="highs")
        if r.status != 0: errs.append("oracle relaxation LP status %d" % r.status)
        elif abs(r.fun - pobj) > 1e-6 * (1 + abs(r.fun)): errs.append("relaxation objective %.9g != oracle optimum %.9g" % (pobj, r.fun))
    return st, len(rows), errs

if __name__ == "__main__":
    try:
        st, k, errs = check(sys.argv[1], sys.argv[2])
    except Exception as ex:  # malformed evidence or an oracle that could not run is a failure, never a pass
        print("status=? FAIL\n  - %s: %s" % (type(ex).__name__, ex)); sys.exit(1)
    print("status=%s rows=%d %s" % (st, k, "PASS" if not errs else "FAIL"))
    for x in errs: print("  -", x)
    sys.exit(1 if errs else 0)
