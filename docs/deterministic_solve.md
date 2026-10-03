# Deterministic work limit for LP solves (opt-in)

`taral MODEL.mps --work-limit N` bounds an LP solve by a cumulative simplex iteration budget instead of
wall-clock time. Without the flag nothing changes: the primal simplex still gets 20 percent of the wall-clock
budget before the dual takes over (`kPrimalShare` in `src/kkt_gate.cpp`), so on a slower or busier host the
route can differ.

## What the flag does

- One counter, `work_used`, counts iterations of the primal and dual simplex engines cumulatively: the primal
  attempt, the dual, the dual's primal clean-up and the equilibrated retry all draw from the same budget.
- Routing uses work, not time: the primal gets `max(1, 20% of N)` iterations; if it stops on that cap and budget
  remains, the dual gets the rest. The wall-clock check for the equilibrated retry is skipped.
- With no `--time-limit` the wall limit is effectively off (1e9 s). If you pass `--time-limit` and it fires, the
  outcome depends on the clock and the determinism claim below does not apply.
- Stopping at the cap returns `iteration_limit` (exit 4) with no objective claim. Status `optimal` still requires
  the existing KKT gate.
- Output: stdout `work_used ... result_hash ...` and JSON fields `work_limit`, `work_used` (clamped to N),
  `result_hash`.
- Scope: continuous linear LP with `--method simplex` or `dual`. `--method ipm`, quadratic objectives, MILP
  (`--node-limit` remains the MILP limit) and `--explain-infeasible` are rejected.

## `result_hash`

FNV-1a 64-bit over canonical text of: status, iteration count, objective (`%.17g`) and every `x` value
(`%.17g`). It excludes wall time and the message (the message carries timings). It is a bit-exact fingerprint,
not a tolerance comparison and not cryptographic.

## What "deterministic" means here (environment and serialization limits)

Claimed: for the same binary, the same input bytes (the MPS text), the same `--work-limit` and the same
method/fallback flags, runs give the same status, iteration count, `work_used` and `result_hash`, regardless of
`--time-limit` and of concurrent machine load. All random perturbations in the engines use fixed seeds and the
engines are single threaded.

Not claimed: identical results across different builds, compilers, optimization flags, CPUs or libm versions.
Floating-point contraction (FMA), instruction selection and library math can change the last bits, and with
pivoting rules the iteration path can then change. Also not claimed: identical hashes after the input is
re-serialized differently (row/column order, number formatting), because those change the parsed model, and
equality of objective values across builds. The hash covers `%.17g` text of doubles produced on one
platform; it is not a portable exchange format.

## Measured (this branch)

`benchmarks/deterministic_solve_tests.py` (all pass): repeated runs, with `--time-limit` unset/5/600, give
identical keys on 6 refinery LP fixtures; identical under 4 busy threads on 3 fixtures; 21 stop-at-cap cases
(caps 1 to 200, 3 fixtures) are repeatable and only ever `iteration_limit` or `optimal`; the hash is the same for
N=1e5 and N=1e7 once the solve completes; flag off reproduces the baseline binary's status, iterations and
objective on the 6 fixtures and writes no work fields; usage and scope errors exit 2 or report `unsupported`.

Cross-build comparison (one host, Intel Xeon 2.6 GHz, g++ 11.4, `--work-limit 100000`, 6 fixtures):
`-O3 -march=native` versus `-O2` versus `-O3 -march=x86-64` versus `-O3 -march=native -ffast-math`. `-O2` and
`-O3 -march=x86-64` agree with each other on all 6; `-O3 -march=native` differs from them on all 6 except
`base_lp_t6`, which matches the `-ffast-math` build only; the `-ffast-math` build differs on the other 5. So the
hash is NOT stable across these build settings, which is why the claim above is limited to one binary.

## Not run / not claimed

- No cross-host run (only one machine type available); no ARM or GPU runs.
- No MILP, IPM or QP determinism work; wall-clock limits remain there.
- No overhead measurement of the counter beyond the test runs above; no iteration-count comparison to the
  default path on large instances (MIPLIB/Netlib not run).
- Whether a given N is a "fair" work unit across models is not claimed: an iteration costs more on larger models.
