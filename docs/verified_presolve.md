# Verified presolve (opt-in)

`taral MODEL.mps --presolve [--presolve-log OUT.json]` reduces a linear LP or MILP before solving, solves the
reduced model with the normal engine, expands the answer to the original columns and audits it in original
space before reporting. Default behaviour (no flag) is unchanged. Not supported: quadratic objectives and
`--method ipm` (usage error).

## Reductions (all floating point, relative tolerance 1e-9)

| op | rule |
|---|---|
| `round_int_bounds` | integer column bounds become ceil(lo) / floor(up): strict inward rounding, NO snap in either direction. An interval with no integer ([4e-7,8e-7], [1.0000004,1.0000008]) is infeasible by exact comparison; an interval containing an integer (including width 0 or 8e-7 around 0, 1, -1, 1e3, 1e9) stays feasible. Without `--presolve` the MILP solver keeps its own documented 1e-6 integrality tolerance and may call such a sub-1e-6 interval optimal; that divergence is recorded, not unified |
| `fix_col` | a column with lo == up is substituted: row bounds and the objective constant shift |
| `empty_row` | a row with no active column is dropped if 0 lies in its bounds, otherwise infeasible |
| `singleton_row` | a one-column row becomes a column bound (integer-rounded); conflict means infeasible |
| `redundant_row` | a row whose min/max activity over the current column bounds lies inside its bounds is dropped |
| `fix_empty_col` | a column in no remaining row is fixed at the bound its cost prefers (skipped if that bound is infinite) |

Infeasibility found by presolve is reported as `infeasible` with the reason in the log and JSON message.
It is a floating-point finding re-derivable from the log, not an exact-rational proof.

## Review fixes (second version)

- Wrong-answer fix: integer bounds [1000000000.4, 1000000000.6] used to round to 1e9 (relative slack 1e-9 x 1e9 = 1) and
  report optimal; now rounding is absolute-slack, the empty interval is `infeasible`, and the original-space bound
  audit uses a capped scale (relative up to 1e3, absolute 1e-6 x 1e3 beyond) so a 0.4 violation at 1e9 fails.
  Regressions: empty integer intervals at 1e9, 1e12 and small scale, plus non-empty and integral big intervals,
  each checked against the plain engine and by the independent replay.
- Fully eliminated models now write the original-space objective and point to `--json`.
- An unwritable or failing `--presolve-log` (fopen, ferror or fclose) exits 5 with `status output_error`.
- Deadline discipline: presolve checks the budget at the start of each pass and every 256 rows/columns. On timeout it
  reports `time_limit` (exit 4), does not use the partial reduction, and writes a log marked `timed_out`. Tested with a
  deterministic hook (`TARAL_PRESOLVE_TEST_TIMEOUT_AFTER_OPS`, test-only) and a C++ API test; a real wall-clock overrun
  was not reproduced or measured (the original finding was static code inspection).
- Replay checker: strict schema validation (integer index ranges per op type, no `-1` where a column is needed, finite
  coefficients and values, NaN/null rejected, audit fields finite and self-consistent), optional `.sol` argument to
  recompute the audit objective, oracle read status checked, duplicate MPS coefficients summed like the engine
  (free-format only; fixed-format names with spaces are unsupported and raise).

## What is verified, and how

- **Replayable log.** `--presolve-log` writes every reduction with the numbers it used. The independent
  `benchmarks/presolve_replay_check.py` reads the ORIGINAL MPS with highspy (not the engine's reader), re-derives
  each reduction from its own state, and requires the kept rows/columns, reduced bounds and objective constant to
  match the log. For an infeasible claim the replayed state must show a violated condition.
- **Original-space audit.** Every point is expanded and checked on the original model by the engine (rows,
  bounds, integrality, recomputed objective incl. the objective constant, vs the solver's reported objective;
  tolerance 1e-6 relative for rows, scaled as above for column bounds, absolute 1e-6 for integrality). A failed audit turns the result into `numerical_failure`, never a reported optimum.
  The same recheck is repeated in the test harness from the `.sol` file.

## Measured (this branch)

`benchmarks/presolve_tests.py` (809 cases at `presolve_tests.py BIN 800`: 5 hand-built, 800 seeded random
LP/MILP with fixed columns, singleton/empty/redundant rows, fractional integer bounds; 4 refinery fixtures):
all pass: plain and presolved status match the oracle, objectives match, every log replays, 464 reported
points audited in original space, 6 corrupted logs rejected (changed fix value, dropped op, extra kept row,
changed constant, changed singleton bound, bogus redundant-row op). Reductions exercised: fix_col 1541,
singleton_row 977, empty_row 387, fix_empty_col 600, round_int_bounds 1255, redundant_row 182.
The four refinery fixtures reduce nothing except `infeasible_supply_lp.mps` (18 fixed columns, 132 -> 114);
the feasible ones are unchanged, so no speed-up is claimed for them.
Existing suites rerun with this binary: `tests/cli_test.sh`, `tests/parse_deadline_test.sh`, and
`benchmarks/milp_tests.py --per-suite 15` (113 ok, 7 undecided at limits, 0 wrong; local run).

## Oracle caveat (test harness)

The oracle is HiGHS through scipy on the model parsed by highspy. HiGHS disagreed with the engine on some random
cases that have integer columns with fractional bounds (for example 0.4) when its own presolve was on, or when the
bounds were passed unrounded. The harness therefore rounds integer bounds itself, turns HiGHS presolve off and
sets `mip_rel_gap` 0. The engine's points in those cases satisfied rows, bounds and integrality in the independent
checker; the cause inside HiGHS is not established here. Integer columns in generated models always have an
explicit upper bound because readers differ on the default.

## Not measured / not claimed

- No MIPLIB comparison. The earlier 26-time-limit MIPLIB ledger is historical, not a baseline, and no new
  MIPLIB run was made, so there is no claim about speed or node counts on real MIPLIB instances.
- No dual postsolve (no duals, reduced costs or bases are mapped back); presolve is for primal answers only.
- No dominated-column, doubleton, coefficient-tightening or probing reductions; no exact arithmetic.
- Reductions are not applied inside branch and bound (root only, `src/milp.cpp` unchanged).
- Comparisons to SCIP/PaPILO/HiGHS presolve are attributed to their documentation only, not measured.

## Output after presolve (primal-only) and replay strictness

- Dual postsolve is not implemented. With `--presolve` the JSON reports the original-space primal point `x`, the
  objective, and `row_activity` recomputed in original space. Reduced-model dual/KKT fields (`dual_objective`,
  `primal_res`, `dual_res`, `gap`, `complementarity`, `max_*`, `row_dual`, `reduced_cost`, residual arrays) are
  `null`, with `original_space_primal_only: true` and `certificate_quality: presolve_primal_only`.
- Infeasible/unbounded results from the reduced model: certificates are not mapped back, `certificate_verified`
  is false, `certificate_quality: presolve_reduced_model_only`, and the message says so. They are reduced-model
  evidence, not original-space verified proofs; the original-space audit applies to reported points only.
  (The presolve-infeasible path with a logged reduction is checked by log replay, not by a certificate.)
- MILP JSON written by the MILP path is unchanged by this change (not re-audited for these fields; unrun).
- Replay comparisons of logged values use a tight tolerance (1e-12 relative, same double arithmetic as the
  engine) and fixed values must equal the bound exactly and be integral for integer columns, so a +1 shift at 1e9
  is rejected. Limit: a shift below about 1e-12 x magnitude is not detected.

## Float-derived integer bounds (singleton rows)

A singleton row bound is computed in floating point (row bound minus fixed-column terms, then divided by the
coefficient), so for an INTEGER column its ceil/floor is not exact. Example: 2X=.394, X+.5Y=.697, Y integer in
[1,10] is feasible (X=.197, Y=1), but float substitution gives .697-.19700000000000006 and a strict floor gave Y<=0,
a false "infeasible". The engine now carries a running absolute error bound per row (one DBL_EPSILON-scaled term per
multiply/subtract, plus the fixed column's own error, plus the division rounding unless |a| = 1) and rounds the
DERIVED integer bound outward by that bound only: ceil(l - e), floor(u + e). Outward rounding only relaxes, so
infeasible claims stay sound; the final answer is still audited in original space. Direct input bounds
(`round_int_bounds`, and singleton rows with |a| = 1 and no substitution) have e = 0 and stay strictly ceil/floor, so
the [4e-7,8e-7] cases remain infeasible. The replay checker recomputes the same error bounds with the same
formulas and rejects the old false-infeasible logs (kept as fixtures). Limits: the error bound is a first-order
model of rounding, not interval arithmetic with directed rounding; an integer bound is treated as exact once rounded.
