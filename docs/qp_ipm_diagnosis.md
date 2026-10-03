# QP and IPM engine gap diagnosis

Branch `cloud/qp-ipm-diagnosis`, started from `origin/cloud/adversarial-tests3` (main `4517277` plus the adversarial harness and the QP suite).
Diagnosis first, then two prototypes on this branch. Nothing here weakens a tolerance, and no `known_failures` file was touched.

Setup: engine built from `src/` with `g++ -O3 -march=native -std=c++17 -Wall -Wextra -Wpedantic` on a 4-core cloud container;
reference HiGHS 1.15.1 (`highspy`); caps 60 s per case (LP) and 30 s (QP); seeds are the suites' own (`26000+N`-style seeds appear
only in the new `tools/diag_qp_ipm/` script). Raw per-case files are in `results/qp_ipm_diagnosis/`.

Two facts about the starting point that differ from the hand-off note:

* `--method` and the IPM ledgers are **not** on `cloud/adversarial-tests3`. They live on the sibling branch `cloud/adversarial-tests2`
  (`tools/adv/harness.py --only ... --method ipm`, `tests/repro_methods/`). Every IPM replay below used a read-only checkout of that branch's
  tools against the engine built from this branch. The QP replays use this branch's `tools/advqp/qpharness.py`.
* The committed ledger counts 68 `dual_infeasible` rows on ipm (56 unbounded + 12 big_status); a fresh run of the same suite on this host gives
  **69**. 13 unbounded cases flip between `unbounded` and `dual_infeasible` between the ledger's build and mine (7 one way, 6 the other, the other 2939
  rows are identical). That flip is itself evidence for gap 2 (a knife-edge status test, below). The QP ledger flips 5 cases the same way
  (`bounds_huge-82057/82058`, `ill_cond-87042/87058` pass in the ledger and fail here, `ill_cond-87054` the other way); on this host
  those are deterministic over repeated runs.

---

## Gap 1: `indef_convex_feasible`, 40 no-answer cases

**Cases:** `indef_convex_feasible-76000` .. `-76039` (all 40). Replay:
`python3 tools/advqp/qpharness.py --engine out/taral --out DIR --cats indef_convex_feasible --jobs 2` (all 40 return `nonconvex`);
per-case table: `python3 tools/diag_qp_ipm/indef_convex_highs.py --engine out/taral`.

**Root cause (engine).** `src/ipm.cpp`, `q_convex()` (called from `ipm_solve` right after the light presolve) runs an LDL^T inertia test on
`Q + 1e-7 max|Q| I` and returns `IpmStatus::Nonconvex` ("Q is not positive semidefinite in the minimisation sense") if any pivot is negative.
It tests the *whole* Q. These models are convex only on the feasible set: `Q` is indefinite (lam_min about -3 in all 40), but `Z'QZ`
(`Z` = null space of the equality rows) is positive definite (lam_min 1.0 .. 10.0). The engine therefore refuses a model that has a unique,
well-defined optimum. This is a deliberate rule ("needs Q PSD"), not a numerical failure, and no solve is attempted.

**HiGHS does answer 20 of the 40, and only because of a cheap check.** HiGHS returns `Optimal` for 20 cases, objective equal to the constructed optimum
to 1e-14 relative (presolve on and off identical), and `Not Set` for the other 20 with the log
`ERROR: Hessian has 10 diagonal entries in [-3, 0) so is not positive semidefinite for minimization / Cannot solve non-convex QP problems with HiGHS`.
Across all 40 the rule is exact: HiGHS refuses iff `Q` has a negative diagonal entry (assertion checked in `diag_qp_ipm`: 40 of 40).
The 20 it answers have an indefinite `Q` whose diagonal happens to be nonnegative, and its active-set solver works in the null space, which is
why it gets the right answer. So HiGHS' 20 are correct but accidental; its other 20 are the same refusal as the engine's. The harness
reference column ("optimal" 20 / "not_set" 20 in `tests3_qp_full.csv`) is exactly this split.

**What is NOT wrong:** no answer is wrong, and the harness does not need a change for these 40 (it already treats them as `feas_convex`).

---

## Gap 2: ipm `dual_infeasible` where HiGHS says unbounded (69 cases)

**Cases:** the 69 are listed by name with before/after status in `results/qp_ipm_diagnosis/ipm_status_gate_changed_rows.csv`
(12 `big_status-63001, 63003 .. 63023` odd ids, 57 `unbounded-36xxx`). Replay (tests2 tools):
`python3 <tests2>/tools/adv/harness.py --engine out/taral --out DIR --only unbounded-36001,big_status-63001 --method ipm --jobs 1`.
Smallest reproducer (`tests/repro_methods/ipm_unbounded_dual_infeasible.mps` on tests2): `MAX 3x s.t. 2x >= 1, x >= 0`, one row, one column.

**It is a weaker true statement, not a wrong label.** `src/ipm.hpp` documents `DualInfeasible` as "recession direction verified, primal
feasibility not established". Independent evidence that all of them really are feasible and unbounded: the generators build them feasible with a ray,
HiGHS (presolve off) says unbounded for all 69, and a zero-objective HiGHS solve of the same constraints is `Optimal` for every one
(`tools/diag_qp_ipm/feasible_check.py`, 69 of 69). The harness scores the status as `noans` on purpose (tests2
docs, finding F1), so it is not a harness mapping bug either.

**Root cause (engine), two separate things in the status path.**

1. `src/ipm.cpp`, `run_ipm()`, at the recession test: `out.status = pres <= 1e-6 ? Unbounded : DualInfeasible`. `pres` is the in-run
   *internal* residual of the infeasible-start iterate, `max(|rp|,|rl|,|ru|)/(1+bnorm)`. It includes the slack-consistency residuals `rl`/`ru`
   of the interior-point formulation and the early iterates are not primal feasible by design. The ray test fires in 1 to 8 iterations. So the status
   depends on a knife-edge on a quantity that is not the feasibility of the point the engine returns. The cross-host flips above are this knife-edge.
2. In `ipm_solve`, after the run the engine already computes the original-model measures of the returned point (`finish()` -> `ms.pres`), and never
   uses them for this decision. Of the 69 base-build points, **40** are feasible on the original model (relative violation `<= 1e-8`, several exactly 0)
   and would be provably unbounded; **29** (all 12 `big_status`, 17 `unbounded`) are not feasible at the iteration where the ray was found
   (`results/qp_ipm_diagnosis/ipm_dual_infeasible_base_points.csv`: original-model relative violations from 1e-8 to 7; two of the 29 are only just above 1e-8). The ray itself is verified only on
   the *scaled* internal problem (`A d = 0`, `Qd = 0` to 1e-9, cone, `c'd < 0`); it is never re-checked on the original model.

The QP suite shows the same thing: all 25 `partial` cases in `unbounded` are this path.

---

## Gap 3: ipm `numerical_failure` on `near_singular` (49) and `huge_bounds` (5)

Ids: `results/qp_ipm_diagnosis/ipm_numerical_failure_ids.txt` (also 1 `near_singular-37146` wrong objective, not analysed).

**near_singular.** The 49 failures are exactly the *plain-cost* variants (`(seed//3)%2 == 0`) of the prescribed-singular-value (cond 1e8..1e12, 23 of 25)
and Hilbert (25 of 25) families, and 1 of 25 of the nearly-dependent-rows family. The objective-stable variants (`c = B'y`) pass 75 of 75.
Trace (`IPM_DBG=1 out/taral case.mps --method ipm`, `near_singular-37002`): dual residual, complementarity and bound residual go to 1e-13 .. 1e-40 while the
primal residual sticks at 3.4e-7 for 30 iterations, `y` keeps growing, then "no progress in 30 iterations". The Newton systems carry static
regularisation `kDelta = 1e-8` (`IPM_DELTA`) and the refinement is against the unregularised matrix; with cond(B) above 1e4 the 1e-8 perturbation is
bigger than the smallest singular value squared, so the solve cannot reduce the primal residual below about 1e-7. Experiment with the engine's own
knobs on the 55 failing ids (no code change; only these 55 were run, not the other 2897):

| setting | near_singular fixed (of the 49 numerical_failure) | huge_bounds fixed (of 5) |
| --- | --- | --- |
| default (`IPM_DELTA=1e-8`) | 0 | 0 |
| `IPM_DELTA=1e-10` | 1 | 0 |
| `IPM_DELTA=1e-12` | 10 | 1 |
| `IPM_DELTA=1e-14` | 20 (8 PASS, 12 PASS_REF: engine certified, reference off) | 0 |
| `IPM_RHO=1e-10` | 0 | 1 |
| `IPM_REF=100` | 0 | 2 |

So static regularisation is a real lever for `near_singular`, but a global change is unsafe (not measured on the rest of the suite).

**huge_bounds.** Reproducer `tests/repro_methods/ipm_huge_bounds_numfail.mps` (tests2): a zero-cost column with a finite upper bound `1e15` in an `E` row of size 5. The starting point centres
inside the box (z about 2e13), nothing pulls it back, and a row of size O(1) cannot be satisfied below the double ulp of 1e13 (about 2e-3): final original residual 7e-4. The scaled problem
reports converged because the primal-residual normaliser `bnorm` (`ipm.cpp`, `const double bnorm`) is the max of `|b|` and **every finite bound**, here 1e15, so the scaled `pres` is about 1e-18
and the run stops with "converged on the scaled problem but not on the original model" after 15 extra iterations.

---

## Fix candidates, ranked by risk (low to high)

Gap 2 (status path):

| # | candidate | risk | state |
| --- | --- | --- | --- |
| D0 | harness: score a verified-ray `dual_infeasible` as `partial` (the QP harness already does) | none for the engine; changes scoring only | not done |
| D1 | promote to `Unbounded` only when the original-model point passes AND the ray passes on the original model (cone per row/column type, `Qd = 0`, `c'd < 0`) | low: status only, optimal path untouched | **prototyped, commit "IPM: promote DualInfeasible ... "** |
| D2 | when the returned point is infeasible, find a feasible point with a zero-objective solve of the same constraints, require the same checks | low to medium: one extra solve (about 3x on the affected cases, no other case) | **prototyped inside D1** |
| D3 | also export the ray / anchor in the JSON for external verification | low, interface change | not done |

Gap 1:

| # | candidate | risk | state |
| --- | --- | --- | --- |
| C0 | keep refusing and relabel the generator expectation in the harness | none for the engine | not done |
| C1 | equality convexification: if `Q + rho A_E'A_E` passes the inertia test (tolerance tied to the original `max|Q|`) the model is convex on `{A_E x = b_E}`; solve the penalised model, which has identical solutions and multipliers there. Only runs when the old test already refused | medium: changes 45 verdicts in the QP suite (below) | **prototyped, separate commit** |
| C2 | full reduced-Hessian check including active inequalities | high, needs a null-space factorisation | not done |
| C3 | inertia-correcting regularisation (local solve of a nonconvex QP) | rejected: gives local optima while the engine reports global certificates | not done |

Gap 3 (not prototyped; one run each only on the failing ids): N1 adaptive `kDelta` (lower it when `dres` and the gap have converged but `pres` is stuck;
evidence above), medium risk because it touches the optimal path; N2 make the stall/convergence test use the original-model measure or cap the bound
term in `bnorm`, medium; N3 start-point clamp for huge finite boxes, medium; N4 a basis-solve polish for square near-singular systems, high effort.

---

## Prototype 1 (D1 + D2): status gate, `src/ipm.cpp`

Rules implemented (all checks on the original unscaled model; any failure keeps `DualInfeasible`):

* the ray comes from the run (`RunOut::ray`, mapped back with the column scaling `D`); max-norm 1; column cone (0 for a boxed column, `>= 0` with only a lower bound,
  `<= 0` with only an upper bound), row cone (`A d` is 0 for a two-sided row, `>= 0` / `<= 0` for one-sided), `Qd = 0`, `c'd < -1e-6 (1 + |c|inf)`;
  slack for the cone and `Qd` checks is `opt.tol` (1e-8), row slack scaled by `sum |a_ij d_j|`;
* a feasible point: the returned iterate if `ms.pres <= opt.tol`, else the point of a zero-objective solve (`c = 0`, `Q = 0`, so it cannot itself end in
  `DualInfeasible` or `Nonconvex`) that passes the same measure;
* no `Infeasible` is ever produced by this path.

Tolerance note: with slack `1e-9` instead of `opt.tol` the gate demotes 4 cases that were `unbounded` before (`unbounded-36004, 36037, 36061, 36102`) and leaves 5 more on
`dual_infeasible`; their ray violations are 1e-10 .. 6e-9 after un-scaling, i.e. iterate noise. The slack is `opt.tol`, the same as the point check;
it is one constant (`rtol` in `ray_ok`) if the program agent wants it tighter.

Before/after, same host, same build flags (`out/` runs; per-case CSVs in `results/qp_ipm_diagnosis/`):

| suite | before | after | moved rows |
| --- | --- | --- | --- |
| ipm LP+big LP (2952, tests2 harness) | PASS 2828, NOANS 123, FAIL_OBJ 1 | PASS 2897, NOANS 54, FAIL_OBJ 1 | 69, all `dual_infeasible` -> `unbounded`; 0 rows with optimal/infeasible status moved; 0 previously-`unbounded` rows moved |
| QP (1780) | 1636 pass, 25 partial | 1661 pass, 0 partial | 25, all `partial` -> `pass`; 0 optimal/infeasible rows moved |
| QP QUICK (200) | 182 pass, 3 partial | 185 pass, 0 partial | the 3 partials |
| base suite QUICK (330, simplex) | | | 0 rows differ |
| Netlib 93, primal and dual, engine run | 92 optimal + dfl001 time_limit / 93 optimal | same | status, objective and iterations identical on all 93 in both modes, except the wall-clock-limited dfl001 iteration count (41464 vs 41231, `time_limit` in both) (compared at full precision; the committed `netlib93_*.txt` files show the objective to 10 significant digits) |
| 5 Farkas cases (`infeasible-35057`, `redundant_eq-34068/34074/34080/34089`) | numerical_failure | numerical_failure | untouched |

Same-status rows are bit-identical in objective and violations (2883 of 2883 on the ipm suite; 209 directly replayed models compared on iteration count too), and
`tests/ipm_gate/check.sh` (new, 6 models) passes: a model that is infeasible *and* has an improving direction stays `dual_infeasible`, never `unbounded`.
Cost: engine time on the 69 cases 10.5 s -> 32.9 s (the zero-objective solves), 206 s -> 219 s over the whole suite. On the near-singular borderline cases of
the safety test below the gate also turns a questionable `unbounded` (strictly convex model, curvature 5e-9) into `dual_infeasible`.

## Prototype 2 (C1): equality convexification, `src/ipm.cpp` (separate commit, gate it separately)

On `{A_E x = b_E}` the objective equals `(c - rho A_E'b_E)'x + 0.5 x'(Q + rho A_E'A_E)x - rho/2 |b_E|^2`; gradient and multipliers of the equality rows are unchanged at feasible points.
`rho` climbs a decade ladder from `max|Q| / max_i |a_i|^2` (9 rungs), the inertia tolerance stays tied to the original `max|Q|` (otherwise the large `rho` would widen it), and the first passing rung is multiplied by 4.
It is attempted only after the old test refused, and gives up if `A_E'A_E` would exceed 2e6 entries.

QP full run on the tree with both prototypes: 1636 -> 1696 pass; `indef_convex_feasible` 0 -> 40/40, `unbounded` 65 -> 90/90. Rows moved: 70 (25 + 40 + 5), nothing else, same-status objectives identical (0 differences).
**Five `indefinite` cases move from PASS to FAIL_STATUS in the harness**: `indefinite-75012, 75027, 75034, 75060, 75076`. The generator labels that category by `lam_min(Q)` only and ignores the
equality rows; for these five `Z'QZ` is positive definite (1 of them, 75034, is PSD with `lam_min` -1e-15), so the model is convex on its feasible set and the engine's answer is correct: objectives agree with
HiGHS on the null-space-reduced convex problem to 1e-9 (`tools/diag_qp_ipm/nullspace_reference.py`, differences 1e-15 .. 1e-9). These are label errors in the suite, not wrong answers, but they count as regressions under the harness as it stands,
so C1 needs either the harness relabel (C0) or a decision to accept them.
Safety test (`tools/diag_qp_ipm/nonconvex_nullspace_negative_test.py`, 300 random models, seed 26007): all 114 models whose `Z'QZ` has `lam_min/|Q|` below -1e-6 are still refused as `nonconvex`;
all 74 with `lam_min` positive (1e-3, 1) are solved, except 1 `numerical_failure`; the 112 borderline models (delta -1e-6, -1e-8, +1e-8) end as `nonconvex` (38, all at -1e-6), `numerical_failure` (70) or, with the gate,
`dual_infeasible` (4); before C1 they were all `nonconvex`. Without the gate, 4 of them were reported `unbounded`, one of them a strictly convex model. QP `qp_repro`: `qp_dual_infeasible_on_unbounded.mps` and `qp_refuses_indefinite_q_convex_on_feasible_set.mps` FIXED, the other four unchanged.

## Reproduce

```
g++ -O3 -march=native -std=c++17 -o out/taral src/*.cpp
python3 tools/advqp/qpharness.py --engine out/taral --out out/qp --jobs 4 --time-limit 30       # QP suite
tests/ipm_gate/check.sh $PWD/out/taral                                                          # status gate models
python3 tools/diag_qp_ipm/indef_convex_highs.py --engine out/taral                              # gap 1 table
python3 tools/diag_qp_ipm/feasible_check.py DIR_WITH_MPS                                        # gap 2 feasibility
python3 tools/diag_qp_ipm/nonconvex_nullspace_negative_test.py out/taral                        # C1 safety test
# IPM-method LP suite: tools from branch cloud/adversarial-tests2
python3 tools/adv/harness.py --engine out/taral --out out/ipm --kind lp+biglp --method ipm --jobs 4 --time-limit 60
```
