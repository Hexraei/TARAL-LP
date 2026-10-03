# Verified presolve (opt-in)

`taral MODEL.mps --presolve [--presolve-log OUT.json]` reduces a linear LP or MILP before solving, solves the
reduced model with the normal engine, expands the answer to the original columns and audits it in original
space before reporting. Default behaviour (no flag) is unchanged. Not supported: quadratic objectives and
`--method ipm` (usage error).

## Reductions (all floating point, relative tolerance 1e-9)

| op | rule |
|---|---|
| `round_int_bounds` | integer column bounds become ceil(lo) / floor(up); a bound within 1e-6 (ABSOLUTE) of an integer snaps to it. The slack is absolute because a relative slack exceeds a whole integer for |bound| above about 1e9. An empty integer interval is infeasible by exact comparison |
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
