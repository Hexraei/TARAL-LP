"""Scoped correction for five legacy indefinite labels, validated in equality space.

No verdict is overridden. Original-model point, objective and KKT checks remain
in qpharness; only classification and the independent reference use this model.
"""
import numpy as np
from qpcase import QCase

CORRECTED_IDS = frozenset('indefinite-%d' % s for s in (75012, 75027, 75034, 75060, 75076))


def reduced_case(case):
    if case.id not in CORRECTED_IDS:
        return None
    eq = np.isfinite(case.rlo) & (case.rlo == case.rup)
    Ae, be = case.A[eq], case.rlo[eq]
    if not len(Ae):
        return None
    _, s, vt = np.linalg.svd(Ae, full_matrices=True)
    rank = int((s > 1e-12 * max(Ae.shape) * max(1.0, s[0])).sum())
    Z = vt[rank:].T
    xp = np.linalg.lstsq(Ae, be, rcond=None)[0]
    if np.max(np.abs(Ae @ xp - be)) > 1e-10 * max(1.0, np.max(np.abs(be))):
        return None
    if Z.shape[1] == 0:
        return None  # fixed-point cases are outside this five-case correction
    H = Z.T @ case.Q @ Z
    H = (H + H.T) / 2
    sign = -1 if case.maximize else 1
    lam = float(np.linalg.eigvalsh(sign * H)[0])
    scale = max(1.0, float(np.linalg.norm(case.Q, 2)))
    if lam < -1e-13 * scale:
        return None
    # Clamp only roundoff-size negative curvature for HiGHS' PSD check.
    if lam < 0:
        H += sign * (-lam) * np.eye(Z.shape[1])
    Aother = case.A[~eq]
    rows = np.vstack([Z, Aother @ Z])
    reduced = QCase('equality_reduced', case.seed, H, Z.T @ (case.c + case.Q @ xp),
                    rows, np.r_[case.lo - xp, case.rlo[~eq] - Aother @ xp],
                    np.r_[case.up - xp, case.rup[~eq] - Aother @ xp],
                    np.full(Z.shape[1], -np.inf), np.full(Z.shape[1], np.inf),
                    maximize=case.maximize, const=case.obj(xp))
    return reduced, lam, Z.shape[1]
