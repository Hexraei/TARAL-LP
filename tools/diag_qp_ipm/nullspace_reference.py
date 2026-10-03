#!/usr/bin/env python3
"""Independent check of the engine's answers on `indefinite` cases whose Q is convex on the equality rows: substitute
x = xp + Z u (Z = null space of the equality rows), which gives a convex QP in u (Hessian Z'QZ) that HiGHS solves, and
compare its objective with the engine's reported one. Cases 12, 27, 34, 60, 76 of the indefinite category (seeds 75012 ...).
  python3 tools/diag_qp_ipm/nullspace_reference.py RESULTS_JSONL   (run from the repo root)"""
import sys, json, numpy as np, highspy
sys.path.insert(0, "tools/advqp")
import qpcase
f = {json.loads(l)["id"]: json.loads(l) for l in open(sys.argv[1])}  # results.jsonl of a qpharness run with the prototype engine
for i in (12, 27, 34, 60, 76):
    case = qpcase.gen("indefinite", i)
    sgn = -1.0 if case.maximize else 1.0
    Q, c = sgn * case.Q, sgn * case.c
    eq = np.isclose(case.rlo, case.rup)
    Ae, be = case.A[eq], case.rlo[eq]
    Ao, rlo, rup = case.A[~eq], case.rlo[~eq], case.rup[~eq]
    n = case.n
    if len(Ae):
        u_, s, vt = np.linalg.svd(Ae); r = int((s > 1e-9).sum()); Z = vt[r:].T; xp = np.linalg.lstsq(Ae, be, rcond=None)[0]
    else:
        Z = np.eye(n); xp = np.zeros(n)
    k = Z.shape[1]
    H = Z.T @ Q @ Z; H = (H + H.T) / 2
    w = np.linalg.eigvalsh(H)[0]
    if w < 0: H = H + (-w + 1e-12) * np.eye(k)   # only the 1e-15 round-off case
    g = Z.T @ (c + Q @ xp)
    rows = np.vstack([Z, Ao @ Z]); rl = np.concatenate([case.lo - xp, rlo - Ao @ xp]); ru = np.concatenate([case.up - xp, rup - Ao @ xp])
    inf = highspy.kHighsInf
    rl = np.where(np.isinf(rl), -inf, rl); ru = np.where(np.isinf(ru), inf, ru)
    lp = highspy.HighsLp(); lp.num_col_, lp.num_row_ = k, rows.shape[0]
    lp.col_cost_ = g; lp.col_lower_ = np.full(k, -inf); lp.col_upper_ = np.full(k, inf)
    lp.row_lower_, lp.row_upper_ = rl, ru
    lp.a_matrix_.format_ = highspy.MatrixFormat.kRowwise
    st, idx, val = [0], [], []
    for r_ in rows:
        nz = np.flatnonzero(r_); idx += nz.tolist(); val += r_[nz].tolist(); st.append(len(idx))
    lp.a_matrix_.start_ = np.array(st, dtype=np.int32); lp.a_matrix_.index_ = np.array(idx, dtype=np.int32); lp.a_matrix_.value_ = np.array(val)
    hs = highspy.HighsHessian(); hs.dim_ = k; hs.format_ = highspy.HessianFormat.kTriangular
    hst, hi, hv = [0], [], []
    for j in range(k):
        for i2 in range(j, k):
            if H[i2, j] != 0: hi.append(i2); hv.append(float(H[i2, j]))
        hst.append(len(hi))
    hs.start_ = np.array(hst, dtype=np.int32); hs.index_ = np.array(hi, dtype=np.int32); hs.value_ = np.array(hv)
    m = highspy.HighsModel(); m.lp_, m.hessian_ = lp, hs
    h = highspy.Highs(); h.setOptionValue("output_flag", False); h.passModel(m); h.run()
    u = np.array(h.getSolution().col_value); x = xp + Z @ u
    ref = case.obj(x)
    ours = f[case.id]["ours_obj"]
    print(case.id, h.modelStatusToString(h.getModelStatus()), "ref obj %.10g" % ref, "engine exact obj", ours, "diff", None if ours is None else abs(ours - ref))
