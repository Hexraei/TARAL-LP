# Per-commit correctness smoke

`bash ci/smoke.sh` builds the current source with the measurement flags and requires zero diagnostics. Python dependencies are test-only, pinned in requirements.txt. CI timings are not benchmark figures.

Coverage:
- Ten SHA-256-verified Netlib models, checked in primal and dual modes against SciPy/HiGHS and an independent original-file row/bound/objective checker. All twenty runs must pass.
- Twelve hand-derived MPS semantics regressions, including ranges, objective constants, name collisions, integers and quadratic input.
- Forty seeded MILP cases and twenty-four seeded QP cases, including infeasible/unbounded models, weak relaxations, limits, general integers, rank deficiency, ill-conditioning and nonconvex rejection. Wrong claims fail. Honest undecided/reference-failed cases are printed separately, not claimed as solved.

This is a small merged-suite adversarial smoke, not the separate cloud adversarial campaign. Add that runner after its branch is reviewed and merged. No full corpus gate, GPU job, schedule, deployment or uploaded artifacts is included.

Workflow runs on pushes, pull requests and manual dispatch on a standard Linux runner. One job, eight-minute hard cap; newer activity on the same branch cancels stale runs. This is a CI signal, not enforced branch protection. Keep the existing account Actions $0 budget/stop-usage setting so exhausted included minutes cannot become paid usage.
