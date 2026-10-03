# Adjudication: QPLIB 8991, engine versus HiGHS reference

Adjudication: the locally reproduced engine point is numerically at the unique global minimum, by positive definiteness and stationarity, and is better than the recorded reference. The raw widening CSV label `WRONG` is retained unchanged; this diagnostic explains the objective discrepancy rather than treating it as an established engine failure.

Problem (QPLIB 8991, converted from the QPLIB .lp file to MPS with highspy as in the widening package): 14,400 continuous variables,
0 constraint rows, box bounds [-0.495867769, 0.495867769], minimise, 71,520 nonzeros in the symmetrised Q.

| Quantity | Value |
|---|---|
| Q positive definite | Gershgorin lower bound 2.0; lam_min by eigsh (tol 1e-6) 2.0006740694 |
| Bound violation of engine point | 0.0 |
| Variables strictly inside their bounds | 14,400 of 14,400 |
| Projected gradient, inf norm, at engine point | 1.12e-12 (largest |c| is 6.83e-4) |
| Engine objective, recomputed from Q, c | -0.0016678673779 |
| 0.5 c'x (equals f at a stationary point of a quadratic) | -0.0016678673776 |
| Reference (HiGHS, harness row) objective | -0.0016618993 |
| Reference minus engine | +5.97e-6 (reference is higher, so worse for a minimisation) |
| Harness rule | WRONG if |engine - ref| > 1e-6 (1 + |ref|) = 1.0017e-6 |

Reasoning: for a positive definite quadratic on a box, a feasible point whose projected gradient is zero is the unique global
minimiser. All 14,400 variables are interior, so the bound constraints are inactive and the point solves Qx = -c. The engine point
satisfies this to 1.1e-12. The reference objective is above the engine objective, so the reference stopped short of the optimum.
The harness flagged the row because the objective is small (1.7e-3) and the rule uses an absolute 1e-6 scale.

Limits of this check: the point and the reference value were reproduced locally on the main-58f77af engine (fresh build) and
compared with the widening-run CSV row (engine -0.0016678674, same to the printed digits); HiGHS was not rerun here (it needs more
memory than the local 2 GB host), so the reference value is taken from the harness row. Do not describe the reference as "wrong",
only as "stopped short of the optimum by 6e-6 absolute".

Reproduce: `python3 adjudicate_convex_qp_box.py 8991.mps engine.sol -0.0016618993` (numpy, scipy, highspy as an MPS reader only).

## Reproduction provenance

The local MPS SHA-256 is `bc9952359da4fc93042f796358d3bc123187d4fcede3e17a3836fad43b330889`, different from the widening-run MPS hash. The original QPLIB LP download hash reported for this re-conversion is `306087d8b6e9a7c19f965112e626e5cbf8ef3af1762dd76162b7138d16c75c5d`, matching the CSV download hash. Original run MPS/point bytes were not retained; mathematical input equivalence is reported from the source hash, not independently proven by a comparison of both serialized models. The local point SHA-256 is `f42a853d694d3c7fac371c209fbb146638ce9c3a921c68beb9d9c87be53fce89`.

The diagnostic was rerun on the attached local model and complete finite 14,400-column point. It confirmed minimization, zero offset, triangular Hessian, Gershgorin lower bound 2.0, all variables interior and projected gradient 1.1207e-12. Since the curvature lower bound is 2, the unconstrained quadratic objective gap at this interior point is bounded numerically by `||gradient||_2^2 / 4`, at most about 4.6e-21 using the recorded infinity norm. This is a floating-point numerical check, not an exact-arithmetic proof file. HiGHS was not rerun.

Run from this folder: `python3 adjudicate_convex_qp_box.py 8991_local.mps 8991_local.sol -0.0016618993072833967`. The script's missing-name default is unsafe for arbitrary point files; the supplied point was separately checked complete with no duplicate names.
