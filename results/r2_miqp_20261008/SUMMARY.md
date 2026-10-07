# R2 convex-MIQP prototype - results package (Oct 8, 2026)

## What was built
Branch-and-bound over QP relaxations (diagonal convex Q, minimisation), built on the
c4800e9 IPM (`ipm_solve`) in a prototype copy of the tree. Scope limits (stated, not hidden):
diagonal convex quadratic only, min sense, no presolve, no warm starts, no cutting planes.

## Corpus and reference
- 1000 random separable convex MIQPs (3-7 integer vars mostly binary, 2-6 continuous,
  3-8 two-sided linear rows), seeds 26119+idx, feasible by construction.
- 32 structured cases across 3 families (mean-variance portfolio with cardinality,
  feature selection with big-M linking, unit commitment with min/max generation),
  seeds 26119000+fam*1000+sub, feasible by construction.
- Reference: exhaustive integer-assignment enumeration; each assignment's continuous
  restriction solved by SLSQP (ftol 1e-12), infeasible assignments skipped. Agreement
  threshold: relative objective difference <= 1e-6 (unchanged from corpus start).

## Results (post-fix, this build)
- Random corpus: 1000/1000 agree. Max rel diff 5.35e-07, mean 5.6e-09.
- Structured corpus: 32/32 agree (11 portfolio, 11 feature-select, 10 unit-commit). Max rel diff 1.06e-08.
- Ledger files: r2_ledger.csv (random), r2_ledger_structured.csv (structured),
  r2_ledger_v1_prepolish.csv (pre-fix provenance, see below).

## Bug found by the corpus and fixed (provenance kept)
Case 409 (seed 26528) exposed a real prototype defect: a near-integer relaxation optimum
(fractionality <= kIntTol = 1e-6) was treated as integer-feasible, so the node was not
branched; the rounded incumbent was 1.31e-06 worse than the true optimum while the
solver's own best_bound (-0.722796741226) sat BELOW its incumbent (-0.722795444366),
contradicting the "optimal" claim. Verified against a tight multi-start reference
(true optimum -0.722796756490 at assignment (0,0,1,1,1), max constraint violation 2.5e-14).
Fix (miqp.cpp, one change): on near-integer nodes, polish by fixing integer columns at
their rounded values and resolving the restricted continuous QP; offer that exact point
(fallback: rounded point). Post-fix case 409: objective -0.722796738381, bound
-0.722796741226, gap 2.8e-09, agrees with reference at 1.8e-08. Full 1000-case corpus
re-run from scratch after the fix: 1000/1000 agree (pre-fix ledger preserved).
Note for the engine lane: src/milp.cpp has a rounding-heuristic-at-every-node with the
same kIntTol = 1e-6 convention; whether a sub-1e-6 analogue of this pattern exists there
is worth a check by the engine agent (I did not touch src/).

## Suggested README replacement for R2
- [x] **R2 - Future solver classes:** convex mixed-integer quadratic programming (MIQP)
  prototype: branch-and-bound over interior-point QP relaxations (diagonal convex Q,
  min sense). Verified 1000/1000 objective agreement against exhaustive-enumeration
  reference on seeded random convex MIQPs and 32/32 on structured portfolio /
  feature-selection / unit-commitment cases (agreement <= 1e-6 relative; measured max
  5.4e-07). Limits: diagonal convex quadratics only; no presolve, warm starts, or cuts;
  nonlinear and mixed-integer nonlinear programming not started.

## Files
- r2harness.py (random-corpus generator + runner + reference), gen_structured.py (3 families)
- Ledgers and generators are archived here; the supplied package did NOT contain the cases/ directory or individual prototype JSON outputs.
- Original measured prototype: c4800e9 copy + miqp.cpp/hpp with a local runner. The original source pair is preserved at ../../experiments/r2_miqp/original/.

## Branch review addendum

The original 1032-case measurements above are NOT a rerun of the branch-review patch.
Review found a separate false-optimal status path: an unresolved relaxation could be
dropped, then a drained queue plus any incumbent reported optimal. The experimental
branch copy now retains inherited bounds for unresolved nodes and reports numerical_failure
or time_limit, even with an incumbent. Failed near-integer polishing is also unresolved.
Eight stubbed-solver status regressions pass; the original copy fails the drained-queue
regression. A standalone build against main 2652581 rechecked random cases 0..9 and
the original near-integer case409: 11/11 agree. Its ledger is branch_smoke_ledger.csv.
No complete 1032-case rerun was performed for this review patch.
The harness header mentions HiGHS QP, but its actual continuous-reference code uses
SciPy SLSQP; the summary and README correctly name exhaustive-enumeration/SLSQP.
This is an experimental numerical prototype, not a certified solver or production CLI.
