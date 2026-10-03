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
import json, sys
import numpy as np, highspy, scipy.sparse as sp
from scipy.optimize import linprog

TOL = 1e-7

def load(path):
    h = highspy.Highs(); h.setOptionValue("output_flag", False)
    h.readModel(path)
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
    e = json.load(open(jpath))
    errs = []
    st = e["status"]
    rows = [r["index"] for r in e["rows"]]
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
            if w is None: errs.append("missing witness for %s" % e["rows"][pos]["name"]); continue
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
        x = np.array(rx["x"]); m, n = A.shape
        lo = np.array([r["lower_relaxed_by"] for r in rx["rows"]]); idx = [r["index"] for r in rx["rows"]]
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
    st, k, errs = check(sys.argv[1], sys.argv[2])
    print("status=%s rows=%d %s" % (st, k, "PASS" if not errs else "FAIL"))
    for x in errs: print("  -", x)
    sys.exit(1 if errs else 0)
