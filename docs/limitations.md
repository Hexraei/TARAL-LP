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
   The post-merge `58f77af` gate measured 93/93 passes in each mode, strict 92/93, with GREENBEA the only nonstrict case ([ledgers](../results/netlib_postmerge_58f77af/)); engine cap 60s/case, pilots counted with HiGHS references capped at 300s, while default `reproduction/reproduce.sh` excludes `pilot.we`/`pilot4` from passes (ceiling 91/93).
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

This selection includes models with 250,997, 1,003,001 and 1,009,306 variables. The widening run below measures this selection; size alone does not establish a timeout or a successful solve. Out-of-class counts describe the library, not engine failures.

## QPLIB measurements

The [19-case widening CSV](../results/qplib_widening_20261003/ledger.csv) records five objective matches and six feasible engine-reported optima without an optimal HiGHS reference. These six are not independently certified optimal solves. Four engine time limits and two numerical failures coincide with reference time limits or solve errors. [Protocol and provenance](../results/qplib_widening_20261003/).

- **8991, reference-gap adjudication:** the raw CSV labels the 5.97e-6 objective discrepancy `WRONG`, but [local reproduction](../results/qplib_widening_20261003/adjudication_8991.md) verifies positive curvature and interior stationarity at the engine point. The recorded reference stopped short by about 6e-6 absolute. This illustrates why a solver-to-reference objective mismatch needs an independent optimality check, not automatic blame. Original run MPS/point bytes were not retained; the local model is re-converted from the reported same-source LP and matches the recorded engine objective.
- **9008, time-cap overrun under investigation:** no JSON after 480.17s under a requested 120s engine cap. The engine investigation identifies a missing deadline check during IPM symbolic factorization setup; a fix is in progress, not verified on main. Other time-limit rows also overrun the requested cap, including 10038 at 140.68s.

Source `58f77af` was verified by matching source-file fingerprints; build and soft caps are recorded with the [provenance](../results/qplib_widening_20261003/provenance.json). Do not extrapolate into a general speed or correctness claim.

## MIPLIB widening measurements

[40 measured cases](../results/miplib_widening_20261003/) on source `58f77af`
record six objective matches, seven infeasibility matches, 26 engine time limits
and one no-JSON anomaly. No row is labeled wrong under this protocol. The
current implementation times out on those 26 cases; HiGHS also reached its
separate 300-second cap on 11 of them. The engine cap was 120 seconds, not an
equal-budget comparison. `neos-1425699` produced no JSON; a crash hypothesis is
under investigation, not confirmed. These are measurements, not 40 solves.
