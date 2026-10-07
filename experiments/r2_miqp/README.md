# Non-shipping diagonal-convex MIQP experiment

This directory is excluded from the production CMake target. The production CLI still
rejects integer quadratic models. Only this standalone runner uses the prototype.
Minimization, diagonal convex Q only; no presolve, cuts, propagation or warm starts.
Relaxation bounds and feasibility use numerical tolerances, not exact certificates.

Build from the repository root:

```sh
experiments/r2_miqp/build.sh /tmp/r2-miqp
/tmp/r2-miqp model.mps --time-limit 60 --json /tmp/result.json
c++ -std=c++17 -Isrc -Iexperiments/r2_miqp experiments/r2_miqp/status_regression.cpp experiments/r2_miqp/miqp.cpp -o /tmp/r2-status-tests
/tmp/r2-status-tests
```

`original/` holds the unchanged measured source. The current copy fixes unresolved-node
status/bound handling and validates relaxation-vector size before access. A failed subtree
or failed near-integer polish cannot be presented as an optimality proof. Incumbents remain
available on numerical failure or time limit. Nonfinite bounds serialize as JSON null.

The status test stubs the QP solver to force paths the measured small corpus did not cover:
numerical failure and time limit after an incumbent, root failure, malformed relaxation,
failed polish, node cap, unsupported maximization, plus a successful control.

[Measured package and review limits](../../results/r2_miqp_20261008/SUMMARY.md).
The archived generators retain their original absolute paths; change PROTO and BASE to
use this runner and a separate scratch directory before rerunning. SciPy is required.
