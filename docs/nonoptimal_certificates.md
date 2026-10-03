# Nonoptimal LP certificates

## Reading the technical terms

LP means linear programming; MILP means mixed-integer linear programming; QP means quadratic programming. IPM is an interior-point method. PDHG is primal-dual hybrid gradient, the approximate GPU method. KKT (Karush-Kuhn-Tucker) checks test feasibility and optimality conditions. A "gate" is the stated validation rule, not an exact-arithmetic proof. "Strict" means the separately stated tighter tolerance; a non-strict pass meets the ordinary rule but not that tighter check. A ledger is a per-case result table. Source hashes, file paths and command flags are retained only so engineers can reproduce a measurement. Historical measurements are not current-source claims.

JSON is the machine-readable result format. Warm starts reuse a previous solution or simplex basis; cold starts do not. fp64 is double-precision floating point; FMA is fused multiply-add.


This document describes infeasibility and unboundedness evidence separately from optimality-condition checks.

For continuous linear models the simplex engine exports original-coordinate proofs in JSON:

- `infeasible`: four nonnegative multiplier arrays, `farkas_row_lower`, `farkas_row_upper`, `farkas_col_lower`, `farkas_col_upper`, in original first-appearance order. They combine the inequalities `-Ax <= -row_lo`, `Ax <= row_up`, `-x <= -col_lo`, `x <= col_up`. A finite-bound weighted contradiction is strictly positive and stationarity cancels. L1 multiplier norm is 1.
- `unbounded`: a feasible original-coordinate `x` anchor and an infinity-normalized `ray`. Finite lower bounds require nonnegative recession activity, finite upper bounds nonpositive activity. The original objective improves along this direction (sense-aware).
- `certificate_verified`, `certificate_residual`, `certificate_margin` describe the internal check. The independent checker ignores these claims and recomputes from the original MPS.

Missing, malformed, nonfinite or failed proof verification returns `numerical_failure`, never a proved nonoptimal status. Time and iteration limits remain unresolved. The dual engine does not export its internal ray: its infeasibility result therefore uses the remaining budget for primal phase 1 to produce a verifiable proof. No time remains means `time_limit`.

This is a floating-point proof contract, not exact rational certification. Stationarity must be <= 1e-12 and any nonzero residual is compensated at its finite minimizing column bound, or a finite bound implied by an original singleton row, computed from the aggregated (row,column) coefficients with exactly one nonzero aggregated column in the row (otherwise rejected); primal/recession violation <= 1e-8; normalized contradiction/improvement must exceed 1e-8. Strict-margin misses stay unresolved. Finite-bound contradictions include inconsistent input bounds. Effective column bounds are checked when the LP API is used inside branch-and-bound. The command-line independent checker is LP-only and rejects integer and quadratic input instead of silently relaxing it. No claim is made about MILP unboundedness certificates or QP rays.

```
g++ -O3 -std=c++17 -Wall -Wextra src/*.cpp -o /tmp/taral
python3 benchmarks/nonoptimal_certificate_tests.py --engine /tmp/taral
g++ -O2 -std=c++17 -Wall -Wextra benchmarks/nonoptimal_certificate_api_tests.cpp src/{mps,simplex,dual,lu,nonoptimal_certificate}.cpp -o /tmp/proof-api
/tmp/proof-api
python3 benchmarks/nonoptimal_certificate_check.py original.mps result.json
```

Synthetic/reference coverage is not a full Netlib gate, CUDA validation or a performance claim. No presolved model, SciPy, or HiGHS is used in the C++ solve/proof path.
