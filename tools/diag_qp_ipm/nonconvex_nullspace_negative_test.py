#!/usr/bin/env python3
"""Safety test for the equality-convexification prototype. 300 random models with free columns and equality rows only,
so convexity is decided by lam_min(Z'QZ) (Z = null space of A), which is set to one of -1, -1e-2, -1e-4, -1e-6, -1e-8, +1e-8,
+1e-3, +1. Prints (class, delta, engine status) counts. A clearly nonconvex model (rel < -1e-6) must be refused as nonconvex.
  python3 tools/diag_qp_ipm/nonconvex_nullspace_negative_test.py ENGINE   (run from the repo root; needs numpy)"""
import sys, os, json, subprocess, collections, numpy as np
sys.path.insert(0, "tools/advqp")
import qpcase
import tempfile
S = tempfile.mkdtemp(); ENG = sys.argv[1]
INF = float("inf")
rs = np.random.RandomState(26000 + 7)  # seed 26007
tab = collections.Counter(); bad = []
for t in range(300):
    n = int(rs.randint(3, 9)); m = int(rs.randint(1, min(4, n - 1) + 1))
    A = rs.randint(-3, 4, (m, n)).astype(float)
    while np.linalg.matrix_rank(A) < m: A = rs.randint(-3, 4, (m, n)).astype(float)
    u, s, vt = np.linalg.svd(A); Z = vt[m:].T; R = vt[:m].T
    delta = [-1.0, -1e-2, -1e-4, -1e-6, -1e-8, 1e-8, 1e-3, 1.0][t % 8]
    k = n - m
    W, _ = np.linalg.qr(rs.randn(k, k)); ev = np.abs(rs.randn(k)) + 1.0; ev[0] = delta
    P = (W * ev) @ W.T
    Q = Z @ P @ Z.T + (R * rs.randn(m)) @ R.T * 3 + 0.7 * (Z @ rs.randn(k, m) @ R.T + R @ rs.randn(m, k) @ Z.T)
    Q = (Q + Q.T) / 2
    x0 = rs.randint(-3, 4, n).astype(float)
    case = qpcase.QCase("neg", 1 + t, Q, rs.randint(-3, 4, n).astype(float), A, A @ x0, A @ x0, np.full(n, -INF), np.full(n, INF))
    p = f"{S}/neg.mps"; open(p, "w").write(qpcase.write_mps(case))
    r = subprocess.run([ENG, p, "--json", f"{S}/neg.json"], capture_output=True, text=True)
    st = json.load(open(f"{S}/neg.json"))["status"]
    lz = float(np.linalg.eigvalsh(Z.T @ Q @ Z)[0]); qm = float(np.abs(Q).max())
    rel = lz / qm
    cls = "nonconvex(rel<-1e-6)" if rel < -1e-6 else ("PSD(rel>1e-6)" if rel > 1e-6 else "border")
    tab[(cls, "delta=%g" % delta, st)] += 1
    if (cls.startswith("nonconvex") and st != "nonconvex") or (cls.startswith("PSD") and st == "nonconvex"):
        bad.append((case.id, cls, st, lz, qm))
for k in sorted(tab): print(k, tab[k])
print('misclassified:', len(bad))
