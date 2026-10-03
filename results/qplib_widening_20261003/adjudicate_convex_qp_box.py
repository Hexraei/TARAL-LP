#!/usr/bin/env python3
"""Adjudicate an engine answer on a convex box-constrained QP (no rows) by certificate, independent of any reference solver.
Usage: python3 adjudicate_convex_qp_box.py model.mps engine.sol [ref_objective]
Needs numpy, scipy, highspy (highspy only as an MPS reader; no solve is run).
Checks: (1) Q symmetric part is positive definite (Gershgorin bound and lam_min by eigsh), (2) point violates no bound,
(3) projected-gradient infinity norm of f at the point, (4) recomputed objective. For a PD quadratic on a box, a point with zero
projected gradient is the unique global minimiser."""
import sys, json, numpy as np, scipy.sparse as sp, highspy
from scipy.sparse.linalg import eigsh
mps, sol = sys.argv[1], sys.argv[2]
ref = float(sys.argv[3]) if len(sys.argv) > 3 else None
h = highspy.Highs(); h.setOptionValue("output_flag", False); h.readModel(mps)
m = h.getModel(); lp = m.lp_; hs = m.hessian_; n = lp.num_col_
assert lp.num_row_ == 0, "rows present: this script handles box-only models"
c = np.array(lp.col_cost_); lo = np.array(lp.col_lower_); up = np.array(lp.col_upper_)
Qm = sp.csc_matrix((hs.value_, hs.index_, hs.start_), shape=(n, n)); Q = (Qm + sp.tril(Qm, -1).T).tocsc()
val = {}
for l in open(sol):
    p = l.split()
    if len(p) == 2: val[p[0]] = float(p[1])
x = np.array([val.get(h.getColName(j)[1], 0.0) for j in range(n)])
f = lambda v: float(c @ v + 0.5 * v @ (Q @ v) + lp.offset_)
rs = np.array(abs(Q).sum(1)).ravel() - abs(Q.diagonal())
g = c + Q @ x
pg = np.where((x <= lo + 1e-12) & (g > 0), 0, np.where((x >= up - 1e-12) & (g < 0), 0, g))
out = dict(n=n, rows=0, nnz_Q_full=int(Q.nnz), gershgorin_min=float((Q.diagonal() - rs).min()),
           lam_min=float(eigsh(Q, k=1, which="SA", return_eigenvectors=False, tol=1e-6)[0]),
           bound_violation=float(max(0, (lo - x).max(), (x - up).max())),
           strictly_interior_vars=int(((x > lo + 1e-12) & (x < up - 1e-12)).sum()),
           projected_gradient_inf=float(abs(pg).max()), c_inf=float(abs(c).max()),
           objective_recomputed=f(x), half_c_dot_x=float(0.5 * c @ x))
if ref is not None:
    out.update(reference_objective=ref, ref_minus_engine=ref - f(x), ref_is_higher_than_engine=bool(ref > f(x)),
               harness_rule_threshold=1e-6 * (1 + abs(ref)), abs_diff=abs(f(x) - ref))
print(json.dumps(out, indent=1))
