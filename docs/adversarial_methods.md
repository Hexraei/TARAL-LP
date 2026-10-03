# Adversarial suite across solver paths (default simplex, dual simplex, interior point) and big models

Everything here is tools and tests; `src/` was not touched. Engine: `src/` at `7d0af4d`, built with `g++ -O3 -march=native -std=c++17`.
Reference: HiGHS 1.15.1 via `highspy` (as in `docs/adversarial.md`). Shared 4-core cloud container, 4 cases at a time, 60 s cap per
engine run and per reference solve; timings are indicative only (one case below changes verdict with the load).

## What was run

`tests/run_adversarial_methods.sh` (one engine build, the same cases and gates for every path):

| path | flag | cases | harness wall |
|---|---|---|---|
| simplex ("before") | none (default primal simplex) | all 3630 base cases + 124 big-model cases = 3754 | 1183 s |
| dual ("after") | `--method dual` | 2840 base LP + 112 big LP = 2952 | 252 s |
| ipm ("after") | `--method ipm` | the same 2952 LP cases | 86 s |

The whole script took 1548 s and **exited 1** (see "Exit status" below). `--method` is only honoured for pure LPs: `src/main.cpp` routes
every MILP to branch and bound whatever the flag says, so the 802 MILP cases are not run on the dual/ipm paths. This was checked, not
assumed: 130 MILP cases (10 per base MILP category) were re-run with `--method dual` and with `--method ipm`; 129 of 130 have the same
verdict and status as the default path in both runs, the 130th (`equality_int-55002`) is `optimal` on all three but scored `PASS` once
and `PASS_REF` twice because HiGHS' presolve answer differs between runs. (Branch and bound itself re-solves nodes with the dual simplex,
`0bad060`, independent of the flag.)

Scoring is unchanged: same status, objective within 1e-6 relative, exact rational feasibility check of the returned point, HiGHS
re-solved with presolve off on any disagreement. The interior-point status `dual_infeasible` is **not** `unbounded` and is scored as no
answer (`noans`), see finding F1.

## New big-model categories

Six categories in `tools/adv/cats/big_*.py` (`BIG = True`: kept out of the base suite and of `tests/seeds.txt`; list with
`python3 tools/adv/gen.py --seeds-big`). Case i of a category uses seed `SEED_BASE + i`, same scheme as the base suite, bases 59000..64000;
the generators use only `random.Random(seed)`. 124 cases, 500 to 3000 rows (`big_assign`: up to 140 rows by 4900 columns).

| category | kind | cases | seeds | what it tests |
|---|---|---|---|---|
| big_sparse | LP | 24 | 59000-59023 | 500-3000 rows, 1.2-2x columns, 4-8 nonzeros per row, every row and column kind, feasible and bounded |
| big_netflow | LP | 24 | 60000-60023 | min-cost flow, 500-3000 nodes, ~4 arcs per node, rank (nodes - 1), cost ties, tight capacities |
| big_staircase | LP | 24 | 61000-61023 | multi-period production planning, inventory chains, shared capacity rows, ranges, free overtime columns |
| big_assign | LP | 16 | 62000-62015 | wide degenerate assignment (N 25-70) and transportation problems |
| big_status | LP | 24 | 63000-63023 | 500-2500 rows made infeasible (Farkas row) or unbounded (recession ray) by construction, alternating |
| big_cover | MILP | 12 | 64000-64011 | set covering, 300-1200 rows, 1.5-3x binary columns |

Engine wall time on the 112 big LP cases (sum over cases / median / worst): simplex 444 s / 0.0-6.2 s per category / 60.2 s;
dual 164 s / 0.0-2.0 s / 31.2 s; ipm 220 s / 0.0-2.9 s / 56.2 s. Two `big_sparse` cases (`59015`, `59021`) exceed the 60 s cap on
the primal simplex path under 4-way load; run alone they finish in 68.6 s and 59.6 s (11840 iterations on the second, about 5 ms per
iteration). That is a slow solve, not a wrong answer, and it is borderline for the cap: with a longer cap they would pass. Dual and ipm solve
both.

## Per-category table (before = simplex, after = dual and ipm)

Cells are `ok / wrong / noans / other` and the pass rate (ok = pass + pass*). wrong = wrong status, objective or an infeasible point;
noans = time or iteration limit, `numerical_failure`, `dual_infeasible`; other = reference unresolved or harness error.

```
category          kind  cases | simplex ok/wrong/noans/other| dual ok/wrong/noans/other | ipm ok/wrong/noans/other  |
---------------------------------------------------------------------------------------------------------------------
degenerate        LP      200 |  200/  0/  0/  0  100.0%  |  200/  0/  0/  0  100.0%  |  200/  0/  0/  0  100.0%  |
cycling           LP      100 |  100/  0/  0/  0  100.0%  |  100/  0/  0/  0  100.0%  |  100/  0/  0/  0  100.0%  |
dup_coef          LP      150 |  150/  0/  0/  0  100.0%  |  150/  0/  0/  0  100.0%  |  150/  0/  0/  0  100.0%  |
dup_empty_zero    LP      200 |  200/  0/  0/  0  100.0%  |  200/  0/  0/  0  100.0%  |  200/  0/  0/  0  100.0%  |
free_fixed        LP      150 |  150/  0/  0/  0  100.0%  |  150/  0/  0/  0  100.0%  |  150/  0/  0/  0  100.0%  |
huge_bounds       LP      150 |  127/ 19/  4/  0   84.7%  |  124/ 19/  7/  0   82.7%  |  145/  0/  5/  0   96.7%  |
negzero           LP       60 |   60/  0/  0/  0  100.0%  |   60/  0/  0/  0  100.0%  |   60/  0/  0/  0  100.0%  |
rank_deficient    LP      150 |  150/  0/  0/  0  100.0%  |  150/  0/  0/  0  100.0%  |  150/  0/  0/  0  100.0%  |
redundant_eq      LP      100 |  100/  0/  0/  0  100.0%  |  100/  0/  0/  0  100.0%  |  100/  0/  0/  0  100.0%  |
infeasible        LP      150 |  150/  0/  0/  0  100.0%  |  150/  0/  0/  0  100.0%  |  150/  0/  0/  0  100.0%  |
unbounded         LP      150 |  150/  0/  0/  0  100.0%  |  150/  0/  0/  0  100.0%  |   93/  0/ 57/  0   62.0%  |
near_singular     LP      150 |  125/ 22/  3/  0   83.3%  |  123/ 25/  2/  0   82.0%  |  100/  1/ 49/  0   66.7%  |
coef_range        LP      200 |  160/ 14/ 26/  0   80.0%  |  187/  4/  9/  0   93.5%  |  200/  0/  0/  0  100.0%  |
ranges            LP      200 |  200/  0/  0/  0  100.0%  |  200/  0/  0/  0  100.0%  |  200/  0/  0/  0  100.0%  |
objconst_sense    LP      120 |  120/  0/  0/  0  100.0%  |  120/  0/  0/  0  100.0%  |  120/  0/  0/  0  100.0%  |
neg_lower         LP      100 |  100/  0/  0/  0  100.0%  |  100/  0/  0/  0  100.0%  |  100/  0/  0/  0  100.0%  |
bound_types       LP      100 |  100/  0/  0/  0  100.0%  |  100/  0/  0/  0  100.0%  |  100/  0/  0/  0  100.0%  |
fuzz_mixed        LP      250 |  250/  0/  0/  0  100.0%  |  250/  0/  0/  0  100.0%  |  250/  0/  0/  0  100.0%  |
tolerance_edge    LP      100 |  100/  0/  0/  0  100.0%  |  100/  0/  0/  0  100.0%  |  100/  0/  0/  0  100.0%  |
large_sparse      LP       60 |   60/  0/  0/  0  100.0%  |   60/  0/  0/  0  100.0%  |   60/  0/  0/  0  100.0%  |
knapsack          MILP     80 |   80/  0/  0/  0  100.0%  | not run (method n/a)      | not run (method n/a)      |
general_int       MILP     80 |   77/  0/  0/  3   96.2%  | not run (method n/a)      | not run (method n/a)      |
setcover          MILP     60 |   60/  0/  0/  0  100.0%  | not run (method n/a)      | not run (method n/a)      |
parity            MILP     40 |   40/  0/  0/  0  100.0%  | not run (method n/a)      | not run (method n/a)      |
fixed_charge      MILP     60 |   60/  0/  0/  0  100.0%  | not run (method n/a)      | not run (method n/a)      |
int_free_neg      MILP     60 |   59/  0/  0/  1   98.3%  | not run (method n/a)      | not run (method n/a)      |
markers_mixed     MILP     80 |   80/  0/  0/  0  100.0%  | not run (method n/a)      | not run (method n/a)      |
int_bound_types   MILP     80 |   80/  0/  0/  0  100.0%  | not run (method n/a)      | not run (method n/a)      |
objconst_max      MILP     50 |   50/  0/  0/  0  100.0%  | not run (method n/a)      | not run (method n/a)      |
equality_int      MILP     60 |   56/  0/  0/  4   93.3%  | not run (method n/a)      | not run (method n/a)      |
unbounded_relax   MILP     40 |   40/  0/  0/  0  100.0%  | not run (method n/a)      | not run (method n/a)      |
big_values        MILP     40 |   38/  0/  2/  0   95.0%  | not run (method n/a)      | not run (method n/a)      |
indicator         MILP     60 |   60/  0/  0/  0  100.0%  | not run (method n/a)      | not run (method n/a)      |
big_sparse        LP       24 |   22/  0/  2/  0   91.7%  |   24/  0/  0/  0  100.0%  |   24/  0/  0/  0  100.0%  |
big_netflow       LP       24 |   24/  0/  0/  0  100.0%  |   24/  0/  0/  0  100.0%  |   24/  0/  0/  0  100.0%  |
big_staircase     LP       24 |   24/  0/  0/  0  100.0%  |   24/  0/  0/  0  100.0%  |   24/  0/  0/  0  100.0%  |
big_assign        LP       16 |   16/  0/  0/  0  100.0%  |   16/  0/  0/  0  100.0%  |   16/  0/  0/  0  100.0%  |
big_status        LP       24 |   24/  0/  0/  0  100.0%  |   24/  0/  0/  0  100.0%  |   12/  0/ 12/  0   50.0%  |
big_cover         MILP     12 |   10/  0/  1/  1   83.3%  | not run (method n/a)      | not run (method n/a)      |
---------------------------------------------------------------------------------------------------------------------
TOTAL             LP     2952 | 2862/ 55/ 35/  0   97.0%  | 2886/ 48/ 18/  0   97.8%  | 2828/  1/123/  0   95.8%  |
TOTAL             MILP    802 |  790/  0/  3/  9   98.5%  | not run (method n/a)      | not run (method n/a)      |
```

pass* counts: simplex 27, dual 29, ipm 0. Totals over the base-suite part only are identical to `docs/adversarial.md` for the simplex path
(55 wrong, 33 noans on the 2840 LP cases); the big LP categories add 110 ok and 2 noans (the two `big_sparse` timeouts).

Failure relation of each path to the default path (per case; machine-readable list in `tests/methods_findings.tsv`):

| path | failing | also fails on simplex | **NEW** (passes on simplex) | fixed (fails on simplex, passes here) |
|---|---|---|---|---|
| dual | 66 | 54 | 12 | 36 |
| ipm | 124 | 19 | 105 | 71 |

## Findings (failures that are new on a path, each with a reproducer or a seed)

Reproducers: `tests/repro_methods/` (replayed with each entry's own `--method`; `tools/adv/repro.py --dir tests/repro_methods`).

**F1. ipm never answers `unbounded` (69 cases).** On every unbounded model where the interior-point path finds the recession
direction it reports `dual_infeasible`, message "recession direction", never `unbounded`: 57 of 150 `unbounded` cases (the other 93 get
`unbounded`) and all 12 unbounded `big_status` cases. Smallest model: `max 3x s.t. 2x >= 1, x >= 0`
(`ipm_unbounded_dual_infeasible.mps`). It is plausibly intended naming (a dual-infeasible LP is primal infeasible or unbounded), so this may
be a status-reporting gap rather than a numerical one; the harness cannot tell a caller anything from it, so it is scored as no answer.
Primal-infeasible detection is fine (12/12 infeasible `big_status` cases, 150/150 `infeasible`).

**F2. ipm gives up on ill-conditioned models.** `numerical_failure` ("no progress in 30 iterations") on 49 `near_singular` cases (30 of them
new; the other 19 also fail on simplex) and 5 `huge_bounds` cases, all feasible and bounded, HiGHS optimal (`ipm_huge_bounds_numfail.mps`;
the `near_singular` ones by seed, see F6). One more `near_singular` case, `near_singular-37146`, is a `FAIL_OBJ`: ipm 6.92585,
HiGHS 6.92583, simplex 6.92563, exact rational optimum 29.99996: all three solvers are far from the exact optimum (condition number
about 1e12) and ipm is outside the 1e-6 gate against HiGHS; simplex and dual are scored `pass*` there only because their point is feasible and
better than HiGHS'.

**F3. dual declares feasible models infeasible (7 new cases).** 6 `near_singular` (`37012, 37018, 37059, 37065, 37113, 37133`) and 1
`coef_range` (`38061`); each time HiGHS has an exactly feasible optimal point (checked in exact arithmetic on the original model). On the
default path these pass. Reproducer for the `coef_range` one: `dual_coef_range_false_infeasible.mps` (2 rows, entries 5e9 and 1.2e-3).

**F4. dual: `numerical_failure` on 4 `huge_bounds` cases** (`31036, 31103, 31123, 31139`), final point violating constraints by 0.93 to 5.6e13, HiGHS optimal;
pass on the default path. Reproducer `dual_huge_bounds_numfail.mps` (1 row, 2 columns, bounds written as the 1e20 sentinel).

**F5. dual: stall on `near_singular-37072`**: `time_limit` after 11 dual iterations and about 20 million primal iterations (hand-over from a stalled dual to
primal simplex, `f40f048`); passes on the default path. Seed only (see F6).

**F6. Reproducer limits for the near_singular findings.** The shrunk models of `37012` (dual false infeasible), `37072` (dual stall) and `37002` (ipm
numerical failure) are 7x6, 9x8 and 4x3. `tools/adv/mkrepro.py` refused them: HiGHS' interior point disagrees with HiGHS' simplex on them, and
an exact rational vertex check of HiGHS' simplex point on these shrunk models did **not** confirm exact feasibility (the active rows are inconsistent in
exact arithmetic), so the shrinker may have turned a feasible near-singular model into one that is infeasible by a hair. They are therefore not installed
and not claimed as reproducers; the full cases are (`python3 tools/adv/harness.py --engine ENGINE --out DIR --only near_singular-37012,near_singular-37072 --method dual`,
`--only near_singular-37002 --method ipm`). On the full cases the returned HiGHS point passes the harness gate (1e-6 relative, exact arithmetic).

**F7. Not new, but changed by the path.** dual fixes 36 cases that fail on simplex (28 `coef_range`, 5 `near_singular`, 2 `big_sparse`, 1 `huge_bounds`) and leaves
the other 54 (19 `huge_bounds` sentinel cases that answer `optimal` on an unbounded model, 15 `near_singular` objective errors, 9 `coef_range` time limits, ...).
ipm fixes 71 (40 `coef_range`, 23 `huge_bounds`, 6 `near_singular`, 2 `big_sparse`): `coef_range` 200/200 and `huge_bounds` 145/150, the best of the three paths there.

**F8. Slow primal simplex at scale.** `simplex_big_sparse_timeout.mps` (`big_sparse-59021`, 1500 rows): 59.6 s alone on the default path, about 2 s on dual,
about 3 s on ipm, about 5 s HiGHS. Replayed with a 30 s cap so it does not depend on load.

## Other observations (not caused by this work)

* `tests/repro/big_values_numfail.mps` replays as `FIXED` (the engine returns optimal -8.782e9, 18 nodes) although `docs/adversarial.md` lists it as a
  numerical_failure; the flag in `tests/repro/expected.json` is stale. The two `big_values` suite cases that time out (`big_values` 38/40) are separate.
* The 4 failures of the default path that are not in `tests/known_failures.tsv` (`big_sparse-59015`, `big_sparse-59021` slow solves; `big_cover-64001` engine time limit where HiGHS
  finishes; `big_cover-64005` neither finishes) are why the default-path report exits 1. I did not add them to the known list: they have no failing-answer reproducer
  except F8.
* `big_cover`: 10 of 12 pass, 1 `noans` (engine time limit, HiGHS optimal), 1 `other` (neither finishes). MILP, so the dual/ipm columns do not apply.
* Seeds: I read the seed instruction as "follow the existing 26000 + 1000k + i scheme", so the new bases are 59000..64000. If a different seed scheme was
  meant, tell me.

## Exit status

`tests/run_adversarial_methods.sh` exit status **1** on the full run (1548 s): harness stages 0, 0, 0; default-path report 1 (4 failures not in the known list, above); before/after
compare 1 (12 new failures on dual, 105 on ipm); `tests/repro` replay 0; method reproducers replayed afterwards by hand (all `KNOWN-FAIL`, 0 regressions). The driver was added
to replay `tests/repro_methods` after the full run started, so that stage is in the script but was not part of the 1548 s run. A QUICK=1 run (10 cases per category, 479 s)
had exited 1 for the same reasons. None of the failures is hidden or added to a known-failures list: exit 0 would need every path-specific failure to be fixed or listed with a reproducer.

## Files

* `tests/run_adversarial_methods.sh`: driver (build, three harness runs, report, compare, both reproducer replays).
* `tools/adv/harness.py`: `--method simplex|dual|ipm`, `--kind big|biglp|lp+big|lp+biglp|all+big`, `--only ID,ID`.
* `tools/adv/compare.py`: the before/after table and NEW/shared/fixed classification; exit 1 on any NEW failure.
* `tools/adv/cats/big_*.py`, `tools/adv/gen.py` (`BIG_CATS`, `--seeds-big`), `tools/adv/core.py` (`rand_core(nnz=...)`, default stream unchanged; base seeds file verified identical).
* `tools/adv/shrink.py`, `repro.py`, `mkrepro.py`: `--method` / per-entry `method` and `time_limit`.
* `tests/repro_methods/`, `tests/methods_findings.tsv`.
