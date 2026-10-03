# Scope and limitations

TARAL-LP is an experimental C++17 solver, not a production replacement for established industrial solvers. An implemented method, a completed measurement and a verified solve are different claims.

## Implemented scope

- LP: CPU primal and dual simplex, an interior-point path, and crossover from an approximate point. The default LP path uses a KKT gate and a primal-budget-then-dual fallback. The fallback gives the first route 20 percent of the available time, so host speed can change the route taken.
- MILP: branch and bound with integer-bound propagation and reduced-cost fixing. Limits and unresolved nodes remain unresolved outcomes, not successful solves.
- QP: continuous convex quadratic objectives with linear or bound constraints through the interior-point path. Returned results require independent checks; implementation does not establish full-library coverage.
- GPU: CUDA PDHG for continuous LPs. It reports `near_optimal`, not an exact vertex or exact-answer certificate. Integer or quadratic models are rejected as `unsupported_model` with exit code 2.

Mixed-integer quadratic models and quadratic constraints are not supported. The current parser rejects unsupported sections and bound types. This is not a claim of complete support for SOS, indicator constraints, nonlinear programming or mixed-integer nonlinear programming.

## Measurement limits

1. No general speed advantage is claimed. A `time_limit`, `iteration_limit`, `node_limit` or numerical failure is not a verified solve.
2. Protocols stay separate. The historical Netlib gate uses 93 denominator cases and a 60-second cap. The pinned Kaggle run at `442ca16` found 90/93 primal passes and 91/93 dual passes, with no wrong answers under that gate. These are historical results, not a statement about every later commit. See the retained [primal](../results/cpp_a9e8218_60s/) and [dual](../results/cpp_0bad060_dual_60s/) ledgers.
3. The Netlib gate checks original-model objective agreement, relative row violation and bound violation. Its pass and stricter residual checks are separate. Other benchmark sets have their own caps, references and checks; completed measurements do not imply optimal solves.
4. No universal "0 wrong" claim is made. Adversarial inputs include known objective discrepancies, timeouts, uncertified cases and unreliable references. Retained adversarial ledgers under `results/adversarial_4517277/` describe that source commit, not current-main performance.
5. Certificates stay strict. A proof the checker cannot verify is refused rather than promoted. Floating-point roundoff can prevent certification even for an instance constructed to be infeasible.
6. GPU and exact CPU comparisons must use the same achieved accuracy. In the retained T4 synthetic scaling table, CPU dual simplex is faster than PDHG on every listed transport size through one million variables. Faster approximate packing answers do not establish faster exact solves. See [the scaling table](../benchmarks/results/gpu_scale_T4_table.csv).
7. Industrial readiness is unverified. Synthetic refinery fixtures are not plant data and do not establish operational savings, integration readiness or consistent million-variable performance.
8. Timing and route choice depend on the machine, compiler flags and load. Near-cap outcomes can move between runs. A reference value alone is not a certificate of feasibility or optimality.

## QPLIB scope

The [QPLIB instance table](https://qplib.zib.de/instances.html), read on October 3, 2026, contains 453 instances. Classification below uses the table's continuous-relaxation convexity, variable type and quadratic-constraint count. The groups are disjoint: quadratic constraints are classified first, then discrete variables, then convexity of the remaining continuous models. Column definitions are in the [QPLIB documentation](https://qplib.zib.de/doc.html).

| Classification | Instances | TARAL-LP scope |
| --- | ---: | --- |
| Continuous, convex objective, no quadratic constraints | 19 | In-class candidates; inclusion is not a completed measurement or solve claim |
| Has quadratic constraints | 284 | Outside implemented scope |
| Remaining models with integer or binary variables | 144 | Outside implemented QP scope |
| Remaining continuous models without the convexity checkmark | 6 | Outside implemented convex-QP scope |

In-class IDs: 10034, 10038, 8495, 8500, 8515, 8547, 8559, 8567, 8602, 8616, 8785, 8790, 8792, 8845, 8906, 8938, 8991, 9002, 9008.

This selection includes models with 250,997, 1,003,001 and 1,009,306 variables. Their performance must be measured; size alone does not establish a timeout or a successful solve. Out-of-class counts describe the library, not engine failures.
