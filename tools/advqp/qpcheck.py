"""Independent checks for QP results. Nothing here uses the engine's reader or the engine's own residuals.

primal_check   exact rational violations of rows and bounds (absolute and 1+|bound| relative) and the exact objective
kkt_certificate  optimality certificate computed from the returned point alone: constraints within tau of a bound
               are taken as active, and an LP over their multipliers y finds the smallest stationarity residual
               of  g = c + Qx = A'y + z  with the KKT sign pattern (z absorbs bound columns). Reported both
               absolute (max_j |r_j|) and relative (max_j |r_j| / (1 + |g_j|)); see docstring for the tolerance.
reference      HiGHS built from the case data with passModel (not from the MPS file), in a forked child with a cap
farkas         independent infeasibility check: HiGHS phase-1 LP (sum of violations) over rows and bounds
"""
import math
import multiprocessing as mp
from fractions import Fraction

import numpy as np
import highspy

INF = math.inf
TAU = 1e-6  # a row/bound is "active" when within TAU * (1 + |bound|) of its bound


def to_frac(a):
    return [Fraction(float(v)) for v in a]


def primal_check(case, x):
    """Exact violations of the generator's own model at x (doubles taken as the exact rationals they are)."""
    n, m = case.n, case.m
    xf = [Fraction(float(v)) for v in x]
    Af = [[Fraction(float(v)) for v in row] for row in case.A]
    ra = rr = ba = br = 0.0
    for i in range(m):
        a = sum((Af[i][j] * xf[j] for j in range(n) if Af[i][j] != 0 and xf[j] != 0), Fraction(0))
        lo, up = case.rlo[i], case.rup[i]
        v, b = Fraction(0), 0.0
        if lo > -INF and a < Fraction(float(lo)):
            v, b = Fraction(float(lo)) - a, lo
        elif up < INF and a > Fraction(float(up)):
            v, b = a - Fraction(float(up)), up
        if v:
            fv = float(v)
            ra, rr = max(ra, fv), max(rr, fv / (1.0 + abs(b)))
    for j in range(n):
        lo, up = case.lo[j], case.up[j]
        v, b = Fraction(0), 0.0
        if lo > -INF and xf[j] < Fraction(float(lo)):
            v, b = Fraction(float(lo)) - xf[j], lo
        elif up < INF and xf[j] > Fraction(float(up)):
            v, b = xf[j] - Fraction(float(up)), up
        if v:
            fv = float(v)
            ba, br = max(ba, fv), max(br, fv / (1.0 + abs(b)))
    # exact objective  c'x + 0.5 x'Qx + const
    obj = Fraction(float(case.const))
    nz = [j for j in range(n) if xf[j] != 0]
    for j in nz:
        obj += Fraction(float(case.c[j])) * xf[j]
    for i in nz:
        for j in nz:
            q = case.Q[i, j]
            if q != 0:
                obj += Fraction(float(q)) * xf[i] * xf[j] / 2
    return dict(row_abs=ra, row_rel=rr, bnd_abs=ba, bnd_rel=br, objective=float(obj))


def _csc(lp, M):
    """Fill lp.a_matrix_ (column-wise) from a dense array."""
    lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
    starts, idx, val = [0], [], []
    for j in range(M.shape[1]):
        nz = np.flatnonzero(M[:, j])
        idx.extend(nz.tolist()); val.extend(M[nz, j].tolist()); starts.append(len(idx))
    lp.a_matrix_.start_ = np.array(starts, dtype=np.int32)
    lp.a_matrix_.index_ = np.array(idx, dtype=np.int32)
    lp.a_matrix_.value_ = np.array(val, dtype=float)


def _solve_lp(lp):
    """Solve a HiGHS LP, retrying with presolve off and then the interior-point solver when simplex reports an error."""
    for opts in ({}, {"presolve": "off"}, {"solver": "ipm", "run_crossover": "off"}):
        h = highspy.Highs()
        h.setOptionValue("output_flag", False)
        h.setOptionValue("threads", 1)
        for k, v in opts.items():
            h.setOptionValue(k, v)
        h.passModel(lp)
        h.run()
        if h.modelStatusToString(h.getModelStatus()) == "Optimal":
            return h
    return None


def dual_certificate(case, x):
    """Optimality certificate from the point alone: multipliers for ALL finite row sides and column bounds
    (y+ >= 0 on lower sides, y- >= 0 on upper sides, z+/z- on bounds), found by two LPs.

      stationarity   min t  s.t.  | (g - A'(y+ - y-) - (z+ - z-))_j | <= t * (1 + |g_j|)      -> stat_rel (stat_abs: unscaled)
      certified gap  min  sum y+*(a_i x - lo_i)^+ + y-*(up_i - a_i x)^+ + z+*(x_j - lo_j)^+ + z-*(up_j - x_j)^+
                     over the multipliers whose stationarity residual is at most max(2 stat_rel, 1e-12)

    With exact stationarity the sum equals f(x) - D(x, y, z) (Wolfe dual value), so for convex Q it bounds
    f(x) - f* from above: gap_rel = gap / (1 + |f(x)|). Minimisation sense (max problems negated first).
    """
    sgn = -1.0 if case.maximize else 1.0
    Qm, cm = sgn * case.Q, sgn * case.c
    x = np.asarray(x, float)
    n, m = case.n, case.m
    g = cm + Qm @ x
    act = case.A @ x if m else np.zeros(0)
    cols, cost_gap = [], []  # per multiplier variable: dense column over the n stationarity equations
    for i in range(m):
        if case.rlo[i] > -INF:
            cols.append(case.A[i, :]); cost_gap.append(max(0.0, act[i] - case.rlo[i]))
        if case.rup[i] < INF:
            cols.append(-case.A[i, :]); cost_gap.append(max(0.0, case.rup[i] - act[i]))
    for j in range(n):
        if case.lo[j] > -INF:
            e = np.zeros(n); e[j] = 1.0; cols.append(e); cost_gap.append(max(0.0, x[j] - case.lo[j]))
        if case.up[j] < INF:
            e = np.zeros(n); e[j] = -1.0; cols.append(e); cost_gap.append(max(0.0, case.up[j] - x[j]))
    k = len(cols)
    M = np.array(cols).T if k else np.zeros((n, 0))  # g - M mult = r
    res = {}
    for mode in ("rel", "abs"):
        w = 1.0 / (1.0 + np.abs(g)) if mode == "rel" else np.ones(n)
        # rows: w_j (g_j - M_j.mult) <= t  and  >= -t  ->  w_j M_j.mult + t >= w_j g_j ; -w_j M_j.mult + t >= -w_j g_j
        Wm = M * w[:, None]
        A = np.vstack([np.hstack([Wm, np.ones((n, 1))]), np.hstack([-Wm, np.ones((n, 1))])])
        lp = highspy.HighsLp()
        lp.num_col_, lp.num_row_ = k + 1, 2 * n
        lp.col_cost_ = np.array([0.0] * k + [1.0])
        lp.col_lower_ = np.zeros(k + 1)
        lp.col_upper_ = np.full(k + 1, highspy.kHighsInf)
        lp.row_lower_ = np.concatenate([w * g, -w * g])
        lp.row_upper_ = np.full(2 * n, highspy.kHighsInf)
        _csc(lp, A)
        h = _solve_lp(lp)
        if h is None:
            res[mode] = (math.nan, None, None, None)
            continue
        t = max(0.0, h.getSolution().col_value[k])
        res[mode] = (t, lp, A, w)
    stat_rel, lp, A, w = res["rel"]
    gap = math.nan
    cap = None
    if lp is not None and k and not math.isnan(stat_rel):
        cg = np.array(cost_gap)
        lp.col_cost_ = np.concatenate([np.where(cg < 1e-13, 0.0, cg), [0.0]])  # distances below 1e-13 are rounding noise
        for cap in [max(2 * stat_rel, 1e-12)] + [10.0 ** (-9 + 0.5 * q) for q in range(13)]:  # loosened only if infeasible
            lp.col_upper_ = np.concatenate([np.full(k, highspy.kHighsInf), [cap]])
            h = _solve_lp(lp)
            if h is not None:
                gap = max(0.0, float(h.getInfo().objective_function_value))
                break
    elif lp is not None:
        gap = 0.0
    if cap is not None and cap > 1e-9 and not math.isnan(gap):
        stat_rel = max(stat_rel, cap)  # the first LP's 0 is within HiGHS' feasibility tolerance; report what was needed
    return dict(stat_rel=stat_rel, stat_abs=res["abs"][0], gap=gap, gap_cap=cap, g_inf=float(np.max(np.abs(g))) if n else 0.0)


def kkt_certificate(case, x, tau=TAU):
    """dual_certificate (gated) merged with active_residual (informational)."""
    out = dual_certificate(case, x)
    out.update(active_residual(case, x, tau))
    return out


def active_residual(case, x, tau=TAU):
    """Stationarity residual when ONLY constraints within tau of a bound may carry multipliers (an exact-active-set
    certificate). An interior-point answer that sits 1e-4 inside a bound that is active at the optimum fails this
    although its objective and duality gap are fine, so it measures how far the point is from the optimal face,
    not whether the answer is optimal. Reported, not gated. Minimisation sense; max problems are negated first.

    Returns stat_rel / stat_abs (best stationarity residual) and n_act (active constraints).
    """
    sgn = -1.0 if case.maximize else 1.0
    Qm = sgn * case.Q
    cm = sgn * case.c
    x = np.asarray(x, float)
    n, m = case.n, case.m
    g = cm + Qm @ x
    act = case.A @ x if m else np.zeros(0)
    # row multipliers: lower-active y>=0, upper-active y<=0, both free, inactive absent
    ycols, ylo, yup, rowidx = [], [], [], []
    for i in range(m):
        lo, up = case.rlo[i], case.rup[i]
        al = lo > -INF and act[i] - lo <= tau * (1 + abs(lo))
        au = up < INF and up - act[i] <= tau * (1 + abs(up))
        if not (al or au):
            continue
        rowidx.append(i)
        ylo.append(-INF if au else 0.0)  # upper-active allows negative
        yup.append(INF if al else 0.0)   # lower-active allows positive
        if al and au:
            ylo[-1], yup[-1] = -INF, INF
    # column sign condition on r = g - A'y : lower-active r_j >= -t, upper-active r_j <= t, both: none, else |r_j| <= t
    cl = np.array([case.lo[j] > -INF and x[j] - case.lo[j] <= tau * (1 + abs(case.lo[j])) for j in range(n)], bool)
    cu = np.array([case.up[j] < INF and case.up[j] - x[j] <= tau * (1 + abs(case.up[j])) for j in range(n)], bool)
    k = len(rowidx)
    Ak = case.A[rowidx, :] if k else np.zeros((0, n))
    out = {}
    for mode in ("rel", "abs"):
        w = 1.0 / (1.0 + np.abs(g)) if mode == "rel" else np.ones(n)
        h = highspy.Highs()
        h.setOptionValue("output_flag", False)
        h.setOptionValue("threads", 1)
        # variables: y_1..y_k, t ; minimise t ;  for each column j:  w_j (g_j - (A'y)_j)  in the allowed band
        lp = highspy.HighsLp()
        lp.num_col_ = k + 1
        lp.col_cost_ = np.array([0.0] * k + [1.0])
        lp.col_lower_ = np.array(ylo + [0.0])
        lp.col_upper_ = np.array(yup + [INF])
        # row j:  -w_j (A'y)_j  -/+ t  <>  -w_j g_j ;  rewrite  r_j = w_j g_j - w_j (A'y)_j  with  -t <= r_j <= t
        # lower-active column: r_j >= -t only (z_j >= 0 absorbs the positive part);  upper-active: r_j <= t only
        free_both = cl & cu
        # r_j + t >= 0 (unless upper-only) and t - r_j >= 0 (unless lower-only)
        # build as two row blocks
        rows_lo = []  # r_j + t >= 0
        rows_up = []  # t - r_j >= 0
        for j in range(n):
            if free_both[j]:
                continue
            if not (cu[j] and not cl[j]):
                rows_lo.append(j)
            if not (cl[j] and not cu[j]):
                rows_up.append(j)
        nr = len(rows_lo) + len(rows_up)
        lp.num_row_ = nr
        A = np.zeros((nr, k + 1))
        lo_b = np.full(nr, -INF)
        up_b = np.full(nr, INF)
        r = 0
        for j in rows_lo:  # w_j g_j - w_j (Ak' y)_j + t >= 0   ->   -w_j Ak[:,j]' y + t >= -w_j g_j
            A[r, :k] = -w[j] * Ak[:, j]; A[r, k] = 1.0; lo_b[r] = -w[j] * g[j]; r += 1
        for j in rows_up:  # t - (w_j g_j - w_j Ak'y) >= 0   ->  w_j Ak[:,j]' y + t >= w_j g_j
            A[r, :k] = w[j] * Ak[:, j]; A[r, k] = 1.0; lo_b[r] = w[j] * g[j]; r += 1
        lp.row_lower_, lp.row_upper_ = lo_b, up_b
        As = A
        lp.a_matrix_.format_ = highspy.MatrixFormat.kRowwise
        starts, idx, val = [0], [], []
        for rr_ in range(nr):
            nzc = np.flatnonzero(As[rr_])
            idx.extend(nzc.tolist()); val.extend(As[rr_, nzc].tolist()); starts.append(len(idx))
        lp.a_matrix_.start_ = np.array(starts, dtype=np.int32)
        lp.a_matrix_.index_ = np.array(idx, dtype=np.int32)
        lp.a_matrix_.value_ = np.array(val, dtype=float)
        h.passModel(lp)
        h.run()
        st = h.modelStatusToString(h.getModelStatus())
        if st != "Optimal":
            out[mode] = (math.nan, None)
            continue
        sol = list(h.getSolution().col_value)
        out[mode] = (max(0.0, sol[k]), np.array(sol[:k]))
    return dict(act_rel=out["rel"][0], act_abs=out["abs"][0], n_act=k + int(cl.sum() + cu.sum()))


# ------------------------------------------------------------------------------------------------ reference
def _hessian_lp(case, sign):
    """HighsModel for  sign * (c'x + 0.5 x'Qx)  (sign=-1 turns a max problem into a min problem)."""
    n, m = case.n, case.m
    lp = highspy.HighsLp()
    lp.num_col_, lp.num_row_ = n, m
    lp.col_cost_ = sign * case.c
    lp.offset_ = sign * case.const
    big = highspy.kHighsInf
    lp.col_lower_ = np.where(np.isinf(case.lo), -big, case.lo)
    lp.col_upper_ = np.where(np.isinf(case.up), big, case.up)
    lp.row_lower_ = np.where(np.isinf(case.rlo), -big, case.rlo)
    lp.row_upper_ = np.where(np.isinf(case.rup), big, case.rup)
    lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
    starts, idx, val = [0], [], []
    for j in range(n):
        nzr = np.flatnonzero(case.A[:, j]) if m else []
        idx.extend(int(i) for i in nzr); val.extend(float(case.A[i, j]) for i in nzr); starts.append(len(idx))
    lp.a_matrix_.start_ = np.array(starts, dtype=np.int32)
    lp.a_matrix_.index_ = np.array(idx, dtype=np.int32)
    lp.a_matrix_.value_ = np.array(val, dtype=float)
    hs = highspy.HighsHessian()
    hs.dim_ = n
    hs.format_ = highspy.HessianFormat.kTriangular
    Qs = sign * case.Q
    hstart, hidx, hval = [0], [], []
    for j in range(n):
        for i in range(j, n):
            if Qs[i, j] != 0:
                hidx.append(i); hval.append(float(Qs[i, j]))
        hstart.append(len(hidx))
    hs.start_ = np.array(hstart, dtype=np.int32)
    hs.index_ = np.array(hidx, dtype=np.int32)
    hs.value_ = np.array(hval, dtype=float)
    mdl = highspy.HighsModel()
    mdl.lp_, mdl.hessian_ = lp, hs
    return mdl


def _solve(case, presolve, tl):
    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    h.setOptionValue("threads", 1)
    h.setOptionValue("time_limit", float(tl))
    h.setOptionValue("presolve", "on" if presolve else "off")
    sign = -1.0 if case.maximize else 1.0
    st0 = h.passModel(_hessian_lp(case, sign))
    if st0 == highspy.HighsStatus.kError:
        return dict(status="ref_load_error", objective=None, x=None)
    h.run()
    ms = h.getModelStatus()
    name = h.modelStatusToString(ms)
    status = {"Optimal": "optimal", "Infeasible": "infeasible", "Unbounded": "unbounded",
              "Primal infeasible or unbounded": "inf_or_unb", "Time limit reached": "time_limit",
              "Not Set": "not_set"}.get(name, "other:" + name)
    out = dict(status=status, objective=None, x=None)
    if status == "optimal":
        x = np.array(h.getSolution().col_value)
        out["x"] = x.tolist()
        out["objective"] = case.obj(x)  # recomputed from the point with the generator's own data
        out["objective_reported"] = sign * h.getInfo().objective_function_value
    return out


def reference(case, presolve=True, tl=30.0, hard=60.0):
    ctx = mp.get_context("fork")
    parent, child = ctx.Pipe(duplex=False)

    def run():
        try:
            child.send(_solve(case, presolve, tl))
        except Exception as e:  # pragma: no cover
            child.send(dict(status="ref_exception:" + repr(e)[:80], objective=None, x=None))

    p = ctx.Process(target=run)
    p.start()
    child.close()
    if parent.poll(hard):
        try:
            out = parent.recv()
        except EOFError:
            out = dict(status="ref_crash", objective=None, x=None)
    else:
        out = dict(status="ref_hang", objective=None, x=None)
    if p.is_alive():
        p.kill()
    p.join()
    return out


def farkas_phase1(case):
    """Minimum total violation of rows and bounds (LP; Q and c ignored). > 1e-9 certifies infeasibility."""
    n, m = case.n, case.m
    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    h.setOptionValue("threads", 1)
    big = highspy.kHighsInf
    # variables x (n), elastic pairs for rows: sp_i, sn_i >= 0 ; for bounds with lo>up we measure lo-up directly
    lp = highspy.HighsLp()
    lp.num_col_, lp.num_row_ = n + 2 * m, m
    lp.col_cost_ = np.array([0.0] * n + [1.0] * (2 * m))
    lo = np.where(np.isinf(case.lo), -big, case.lo)
    up = np.where(np.isinf(case.up), big, case.up)
    crossed = float(np.sum(np.maximum(0.0, case.lo - case.up)))
    lo2 = np.minimum(lo, up)  # keep the LP loadable; the crossed amount is added back below
    up2 = np.maximum(lo, up)
    lp.col_lower_ = np.concatenate([lo2, np.zeros(2 * m)])
    lp.col_upper_ = np.concatenate([up2, np.full(2 * m, big)])
    lp.row_lower_ = np.where(np.isinf(case.rlo), -big, case.rlo)
    lp.row_upper_ = np.where(np.isinf(case.rup), big, case.rup)
    lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
    starts, idx, val = [0], [], []
    for j in range(n):
        nzr = np.flatnonzero(case.A[:, j]) if m else []
        idx.extend(int(i) for i in nzr); val.extend(float(case.A[i, j]) for i in nzr); starts.append(len(idx))
    for i in range(m):  # row_i + sp_i - sn_i
        idx.append(i); val.append(1.0); starts.append(len(idx))
        idx.append(i); val.append(-1.0); starts.append(len(idx))
    lp.a_matrix_.start_ = np.array(starts, dtype=np.int32)
    lp.a_matrix_.index_ = np.array(idx, dtype=np.int32)
    lp.a_matrix_.value_ = np.array(val, dtype=float)
    h.passModel(lp)
    h.run()
    if h.modelStatusToString(h.getModelStatus()) != "Optimal":
        return math.nan
    return float(h.getInfo().objective_function_value) + crossed


def verify_ray(case):
    """Exact check of the constructed recession direction: A d in the row cone, d in the column cone, Qd=0, c'd<0
    in the objective's own sense. Returns an error string, or '' if the construction is sound."""
    d = [Fraction(float(v)) for v in case.ray]
    n, m = case.n, case.m
    for j in range(n):
        if d[j] > 0 and case.up[j] < INF:
            return "ray leaves the upper bound of column %d" % j
        if d[j] < 0 and case.lo[j] > -INF:
            return "ray leaves the lower bound of column %d" % j
    for i in range(m):
        ad = sum((Fraction(float(case.A[i, j])) * d[j] for j in range(n)), Fraction(0))
        if ad > 0 and case.rup[i] < INF:
            return "ray leaves row %d from above" % i
        if ad < 0 and case.rlo[i] > -INF:
            return "ray leaves row %d from below" % i
    for i in range(n):
        if sum((Fraction(float(case.Q[i, j])) * d[j] for j in range(n)), Fraction(0)) != 0:
            return "Q d != 0"
    cd = sum((Fraction(float(case.c[j])) * d[j] for j in range(n)), Fraction(0))
    if (cd >= 0) if not case.maximize else (cd <= 0):
        return "c'd does not improve the objective"
    return ""
