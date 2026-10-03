# Infeasible-model sweep (row-irreducible explanations)

`benchmarks/infeasible_sweep.py` runs the engine's plain solve and `--explain-infeasible` over a directory of
MPS files plus seeded generated LPs, checks every explanation with the independent checker
(`benchmarks/infeasibility_explanation_check.py`, HiGHS as test oracle only) and optionally records HiGHS' IIS size.
No corpus is stored in the repository: pass `--corpus DIR`. The measurements below used the 29 Netlib infeasible
models supplied as a package (archive SHA256 9e8e1870187a613591a1e41addd2ba6f44162bf9da93402c943017ac8afdf7db, every
file matched the package manifest's `mps_sha256`; not committed here).

## What each procedure guarantees

| procedure | guarantee | NOT guaranteed |
|---|---|---|
| taral `irreducible` | row-irreducible under retained column bounds: the rows plus all column bounds have a verified Farkas certificate (floating point, 1e-8 margin) and removing any one row leaves a stored witness feasible at 1e-7 | minimum cardinality, uniqueness, exact-rational proof |
| taral `reduced_unproven` | rows plus bounds have a verified Farkas certificate; some rows' removal tests were inconclusive (`unproven_rows`) or timed out | irreducibility; the set may contain removable rows |
| HiGHS IIS (`iis_strategy` 6 = FromLp 2 + Irreducible 4, HiGHS 1.15.1) | whatever HiGHS documents for those flags; sizes are reference values only | no claim of irreducibility or minimality is made here from the counts; rows and columns are counted together (a column entry is a bound), so they are not row-for-row comparable with the taral row count |
| planted generated LPs | the planted conflict size k is an UPPER bound on the smallest conflict | a size equal to k proves nothing about minimality |

Nothing here is labelled "smallest". A provably smallest conflict set is a separate, not-yet-built search.

## Measured (this branch, one host; the HiGHS IIS columns come from one run)

| model | engine | explain | rows | unproven | HiGHS IIS rows/cols | explain s | row-proof check | relaxation check |
|---|---|---|---|---|---|---|---|---|
| bgdbg1 | infeasible | irreducible | 5 | 0 | 3/29 | 0.016 | PASS | PASS |
| bgetam | infeasible | reduced_unproven | 29 | 11 | 7/24 | 0.062 | PASS | PASS |
| bgindy | infeasible | irreducible | 3 | 0 | 3/993 | 2.838 | PASS | PASS |
| bgprtr | infeasible | irreducible | 5 | 0 | 5/9 | 0.003 | PASS | PASS |
| box1 | infeasible | irreducible | 114 | 0 | 8/11 | 0.069 | PASS | PASS |
| ceria3d | infeasible | reduced_unproven | 808 | 2 | 73/93 | 37.892 | PASS | PASS |
| chemcom | infeasible | irreducible | 7 | 0 | 7/37 | 0.043 | PASS | PASS |
| cplex1 | infeasible | irreducible | 5 | 0 | 5/6 | 4.071 | PASS | FAIL |
| cplex2 | numerical_failure | not_run_engine_not_infeasible | - | - | 224/221 | - | - | - |
| ex72a | infeasible | irreducible | 152 | 0 | 58/68 | 0.175 | PASS | PASS |
| ex73a | infeasible | irreducible | 103 | 0 | 24/28 | 0.07 | PASS | PASS |
| forest6 | infeasible | irreducible | 57 | 0 | 59/95 | 0.028 | PASS | PASS |
| galenet | infeasible | irreducible | 2 | 0 | 2/4 | 0.002 | PASS | PASS |
| gosh | numerical_failure | not_run_engine_not_infeasible | - | - | 9/8 | - | - | - |
| gran | infeasible | reduced_unproven | 1035 | 1035 | 1/1 | 60.026 | PASS | PASS |
| greenbea | infeasible | irreducible | 1 | 0 | 1/93 | 3.924 | PASS | PASS |
| itest2 | infeasible | irreducible | 3 | 0 | 3/3 | 0.002 | PASS | PASS |
| itest6 | infeasible | irreducible | 2 | 0 | 2/3 | 0.002 | PASS | PASS |
| klein1 | numerical_failure | not_run_engine_not_infeasible | - | - | 50/54 | - | - | - |
| klein2 | numerical_failure | not_run_engine_not_infeasible | - | - | 60/54 | - | - | - |
| klein3 | numerical_failure | not_run_engine_not_infeasible | - | - | 90/88 | - | - | - |
| mondou2 | infeasible | irreducible | 32 | 0 | 26/70 | 0.03 | PASS | PASS |
| pang | infeasible | reduced_unproven | 12 | 2 | 11/35 | 0.022 | PASS | PASS |
| pilot4i | infeasible | irreducible | 1 | 0 | 1/42 | 0.15 | PASS | PASS |
| qual | infeasible | reduced_unproven | 135 | 58 | 79/173 | 0.708 | PASS | FAIL |
| reactor | infeasible | irreducible | 1 | 0 | 1/14 | 0.044 | PASS | PASS |
| refinery | infeasible | irreducible | 51 | 0 | 54/126 | 0.099 | PASS | PASS |
| vol1 | numerical_failure | not_run_engine_not_infeasible | - | - | 84/176 | - | - | - |
| woodinfe | infeasible | irreducible | 1 | 0 | 1/1 | 0.003 | PASS | PASS |

Counts on the 29 Netlib models (engine on main ba03938, `--time-limit 60`):
- engine status `infeasible`: 23. engine status `numerical_failure`: 6 (cplex2, gosh, klein1, klein2, klein3, vol1).
  This is the verdict of the existing engine; these six are NOT explained and nothing is claimed about them. HiGHS
  reports Infeasible for all 29, but the engine's verdict is not certified for the six, and 29 are not claimed certified.
- of the 23 infeasible: 18 explanations `irreducible`, 5 `reduced_unproven` (bgetam, ceria3d, gran, pang, qual).
  `gran` hit the 60 s limit with all 1035 rows unproven.
- the independent checker accepted the row-proof part (rows plus bounds infeasible, per-row witnesses) for all 23.
- the relaxation part (elastic LP optimum compared with a HiGHS simplex optimum at relative 1e-6) failed on two
  models: cplex1 (engine 0.3826807 vs HiGHS simplex 0.377664455; HiGHS interior point gives a third value, 0.39088, so
  the optimum of this badly scaled LP is not independently confirmed) and qual (0.00541016 vs 0.00542051). The engine
  labelled both relaxations `optimal`. That label is a tolerance-level KKT statement, not an independently confirmed
  optimum, on these two models. Not fixed here (solver-core change, off this lane).
- the engine's irreducible sets can be much larger than a HiGHS IIS: box1 114 rows vs 8, ex72a 152 vs 58, ex73a 103 vs 24.
  Row-irreducible is a weak notion of small. (HiGHS counts above are reference sizes only.)

Generated LPs (`--generated 60 --seed 7`, neutral seeds): 60 of 60 explained `irreducible` and accepted by the checker;
the explained set was no larger than the planted k in 60 of 60 and equal to k in 58 (an upper bound, not a minimality
statement). Smaller HiGHS IIS than the engine's set: 0 of 60; larger: 0 of 60 (row counts).

## Checker change in this branch

The independent checker's oracle model was wrong for `greenbea`: the file defines several named BOUNDS sets; the engine
uses only the first set and lets the last entry per side win, HiGHS keeps the first entry and warned. The checker now
normalizes the BOUNDS section the way the engine reads it (tested on small cases). Before the change the checker
could not read greenbea; after it, the row proof is accepted. Only LO/UP/FX/MI/PL/FR are normalized; BV/LI/UI/SC bounds
are passed through unchanged.

## Not run / not claimed

- No provably-minimum-cardinality search (separate design). No cross-host run. MIPLIB infeasible models not run.
- No explanation of the six `numerical_failure` models, by instruction, until the engine verdict itself is settled.
- Exact-rational verification of certificates is not done; certificates are floating point.
- Wall times are from one host and were not repeated.
