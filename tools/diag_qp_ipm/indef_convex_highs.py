#!/usr/bin/env python3
"""indef_convex_feasible diagnosis: for each of the 40 cases, lam_min(Q), lam_min(Z'QZ) (Z = null space of A), number of
negative diagonal entries of Q, HiGHS status (presolve on and off) and, when HiGHS answers, its distance from the
constructed optimum; with --engine, also the engine status. Needs numpy and highspy; run from the repo root.
  python3 tools/diag_qp_ipm/indef_convex_highs.py [--engine out/taral]"""
import json, os, subprocess, sys, tempfile
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "advqp"))
import qpcase, qpcheck

engine = sys.argv[sys.argv.index("--engine") + 1] if "--engine" in sys.argv else None
tmp = tempfile.mkdtemp()
print("id n m lam_min(Q) lam_min(Z'QZ) neg_diag highs_presolve_on highs_presolve_off highs_obj_vs_construction engine")
for i in range(40):
    case = qpcase.gen("indef_convex_feasible", i)
    A, Q = case.A, case.Q
    Z = np.linalg.svd(A)[2][A.shape[0]:].T
    hs = []
    for pres in (True, False):
        r = qpcheck.reference(case, presolve=pres, tl=20)
        hs.append(r)
    nd = int((np.diag(Q) < 0).sum())
    d = "-" if hs[0]["status"] != "optimal" else "%.1e" % (abs(hs[0]["objective"] - case.expect_obj) / max(1, abs(case.expect_obj)))
    est = "-"
    if engine:
        p = os.path.join(tmp, "c.mps"); open(p, "w").write(qpcase.write_mps(case))
        subprocess.run([engine, p, "--json", p + ".json"], capture_output=True)
        est = json.load(open(p + ".json"))["status"]
    print(case.id, case.n, case.m, "%.3f" % np.linalg.eigvalsh(Q)[0], "%.3f" % np.linalg.eigvalsh(Z.T @ Q @ Z)[0], nd,
          hs[0]["status"], hs[1]["status"], d, est)
