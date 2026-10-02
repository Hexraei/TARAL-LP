"""Independent check of a PDHG answer on the ORIGINAL model, sharing no code with gpu/pdhg.cu or src/.

Uses benchmarks/orig_check.py (its own MPS parser) for the primal side: row and bound violation and the
objective, incl. the objective-row RHS constant. Dual side: any row multipliers y give a weak-duality lower
bound for  min c'x  s.t.  rows in [lo,up], x in [l,u]:
    LB(y) = sum_j min_{x_j in [l_j,u_j]} d_j x_j  +  sum_i min_{r_i in [lo_i,up_i]} y_i r_i  + const,   d = c - A'y.
LB is a valid bound for ANY y (either sign convention). A term with an infinite bound and a wrong-sign coefficient
makes LB -inf, which an approximate PDHG dual nearly always does. So three numbers are reported: chk_dual_bound (the
strict weak-duality bound, empty when -inf), chk_dual_obj_loose (those terms dropped: NOT a rigorous bound) and
chk_dual_infeas_rel (the largest dropped term, relative). The sign (+y or -y) with the smaller infeasibility is used.
Models with an OBJSENSE section are not dual-checked (the parser reads minimisation only).
"""
import math, os, sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'benchmarks'))
import orig_check as oc
from pathlib import Path

INF = 1e20


def _row_bounds(model, r):
    k = model['rows'][r]; b = model['rhs'].get(r, 0.); R = model['rng'].get(r)
    if k == 'L': return (b - abs(R) if R is not None else -INF), b
    if k == 'G': return b, (b + abs(R) if R is not None else INF)
    if R is None: return b, b
    return (b, b + R) if R >= 0 else (b + R, b)


def _min_over(coef, lo, up):
    """(min of coef*v over [lo,up], infeasibility): an infinite side with the wrong sign gives -inf; the term is then
    dropped from the loose sum and its size recorded as dual infeasibility."""
    if coef == 0: return 0.0, 0.0
    v = lo if coef > 0 else up
    return (None, abs(coef)) if abs(v) >= INF else (coef * v, 0.0)


def _read_vec(path):
    d = {}
    for line in open(path):
        k, v = line.rsplit(None, 1)
        d[k] = float(v)
    return d


def check_files(mps, sol, dual=None):
    model = oc.parse_mps(Path(mps))
    out = oc.check(model, _read_vec(sol))
    res = dict(chk_objective=out['objective'], chk_row_viol_rel=out['row_violation_rel'], chk_bound_viol=out['bound_violation'],
               chk_dual_bound=None, chk_dual_obj_loose=None, chk_dual_infeas_rel=None, chk_dual_gap_rel=None, chk_dual_sign=None)
    if dual is None or 'OBJSENSE' in Path(mps).read_text():
        return res
    y0 = _read_vec(dual)
    const = -model['rhs'].get(model['obj'], 0.)
    best = None
    for sgn in (1.0, -1.0):
        y = {r: sgn * y0.get(r, 0.) for r in model['order']}
        loose, strict_ok, infeas = const, True, 0.0
        for r in model['order']:
            t, inf = _min_over(y[r], *_row_bounds(model, r))
            if t is None: strict_ok = False; infeas = max(infeas, inf / (1. + abs(y[r])))
            else: loose += t
        for j in model['corder']:
            c = model['cols'][j].get(model['obj'], 0.)
            d = c - sum(a * y[r] for r, a in model['cols'][j].items() if r in y)
            t, inf = _min_over(d, *model['bnd'].get(j, [0., INF]))
            if t is None: strict_ok = False; infeas = max(infeas, inf / (1. + abs(c)))
            else: loose += t
        if best is None or infeas < best[0]: best = (infeas, sgn, loose, strict_ok)
    infeas, sgn, loose, strict_ok = best
    res.update(chk_dual_sign=sgn, chk_dual_infeas_rel=infeas, chk_dual_obj_loose=loose,
               chk_dual_bound=loose if strict_ok else None,
               chk_dual_gap_rel=(out['objective'] - loose) / max(1.0, abs(out['objective'])))
    return res
