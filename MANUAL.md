# TARAL-LP evaluator manual

This manual is for technical evaluators seeing the repository for the first time.
It covers building the solver from source, running it, running the test suites,
finding your way around the tree, and making and checking your own edits. It
describes mechanisms, not results: every measured claim lives in the README and
the documents under `docs/`, each linked to its evidence under `results/`.
Where this manual and `taral --help` disagree, `--help` wins for that source
snapshot.

## 1. Prerequisites

| Need | For |
| --- | --- |
| A C++17 compiler (`g++` is what the recorded measurements used) | building the solver |
| CMake 3.16 or newer for builds; 3.20 or newer for `scripts/dev.sh test` (it uses `ctest --test-dir`) (optional) | the CMake build and test routes |
| Python 3.10+ with `pip` | benchmark and check tooling |
| CUDA toolkit (optional) | the experimental GPU code under `gpu/` only |

The CPU solver itself uses only the C++ standard library. Python packages
(HiGHS, NumPy, SciPy) are used by benchmark and verification tooling, never by
the optimization engine.

## 2. Build from source

CMake route, from the repository root:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

This produces the solver as `build/taral`. A plain `cmake` configure without
`-DCMAKE_BUILD_TYPE` applies no optimization flags; pass `Release` for speed.
Two optional CMake switches (both default OFF): `-DTARAL_WARNINGS=ON` adds
`-Wall -Wextra -Wpedantic`, and `-DTARAL_BUILD_TESTS=ON` builds the C++ unit
tests and registers them with `ctest` (see section 4).

The same routes are wrapped by one developer entry point:
`scripts/dev.sh build` (CMake Release with warnings, `build/taral`),
`scripts/dev.sh pinned` (the exact measurement build used for the committed
ledgers, producing `./taral`), and `scripts/dev.sh smoke` (section 4).

Single-command route (the form used for the recorded interface measurements):

```sh
g++ -O3 -std=c++17 -o taral src/*.cpp
```

The CI smoke builds with stricter flags and requires a warning-free compile:

```sh
g++ -O3 -march=native -std=c++17 -Wall -Wextra -Wpedantic -o taral src/*.cpp
```

There is no install step and there are no third-party build dependencies.

## 3. Run the solver

```sh
./taral MODEL.mps --time-limit 30 --sol /tmp/out.sol --json /tmp/out.json
./taral --help
```

Input is a model file in MPS format. `--help` prints the full option list for
the build you are running; `docs/interface.md` documents the options and the
measured interface scope. Commonly used options:

| Option | Meaning |
| --- | --- |
| `--time-limit S` | Wall-time budget in seconds (default 60); setup can overrun it, it is not a hard deadline |
| `--node-limit N` | Branch-and-bound node limit (default unlimited) |
| `--method simplex\|dual\|ipm` | Continuous-method selection (default primal simplex) |
| `--work-limit N` | Deterministic iteration budget (integer >= 1) for continuous linear LPs under simplex or dual only; cannot be combined with `--presolve` or `--explain-infeasible`, does not apply to IPM, QP or MILP; an explicit `--time-limit` that fires still stops the run |
| `--presolve`, `--presolve-log F` | Verified presolve with a replayable log and original-space audit |
| `--sol FILE`, `--json FILE` | Solution and machine-readable result destinations |

Exit codes: 0 for a definitive status (optimal, infeasible, unbounded), 2 for a
usage error, 3 for a parse error, 4 for a time/iteration/node limit, 5 for
other failures. Two habits matter when evaluating:

- A file existing is not a result. Inspect the status inside the JSON or
  standard output, and treat a time-limited mixed-integer run as carrying an
  unproven incumbent (`has_solution`, `best_bound`, `gap`), not an optimum.
- A definitive status is a floating-point result under the validation rules
  stated in `docs/measured-results.md`, not an exact-arithmetic proof.

## 4. Run the test suites

Ordered from fastest to most involved. All commands run from the repository
root, assume the solver binary at `build/taral` or `./taral`, and need no
network access unless stated.

| Suite | Command | What it covers |
| --- | --- | --- |
| Unit + CLI tests (ctest) | `scripts/dev.sh test` | Builds with `-DTARAL_BUILD_TESTS=ON` and runs `ctest` over the 10 wired tests: the C++ unit tests (`tests/kkt_gate_tests.cpp`, the `benchmarks/*_api_tests.cpp` API tests), the CLI and regression shell scripts, and the IPM gate |
| One-shot build + smoke | `bash scripts/build_and_smoke.sh` | CMake build, then small deterministic benchmark checks in a fresh Python venv; CPU only and no benchmark-corpus downloads, but package installation requires network access or an offline wheel cache |
| Command-line validation | `sh tests/cli_test.sh build/taral` | Help text, valid and invalid options, exit codes |
| Presolve/work-limit rejection | `sh tests/presolve_work_limit_test.sh build/taral` | Confirms `--presolve` combined with `--work-limit` is rejected with a usage error (exit 2) before execution |
| Parse deadline | `sh tests/parse_deadline_test.sh build/taral` | Parser behavior under time limits (the script defaults to `./taral`; pass the binary explicitly after a CMake build) |
| Python benchmark tests | interfaces differ per script; each docstring is the contract | Examples: `python3 benchmarks/deterministic_solve_tests.py build/taral` (positional binary), `python3 benchmarks/milp_tests.py --engine build/taral` (`--engine` flag), `python3 benchmarks/rim_vector_report.py model.mps` (static report; takes model files, no solver) |
| Per-commit CI smoke | `python3 -m pip install -r ci/requirements.txt`, then `bash ci/smoke.sh` | Warning-free build, ten hash-verified Netlib models in primal and dual modes, MPS-semantics regressions, seeded MILP and QP cases; downloads a small pinned corpus on first run |
| Full reproduction | `python3 -m pip install -r reproduction/requirements.txt`, then `bash reproduction/reproduce.sh` | Fetches the pinned Netlib corpus, verifies every file hash, builds the engine, runs the 93-case gate plus TRUSS (TRUSS sits outside the 93-case denominator). `pilot.we` and `pilot4` are both run but excluded from counted passes under the 60 s reference protocol (ceiling 91/93); see `reproduction/README_REPRO.md` |

Python tooling dependencies for the smoke route are pinned in
`benchmarks/requirements-smoke.txt`, and the CI smoke pins its own set in
`ci/requirements.txt`. The reproduction route's `reproduction/requirements.txt`
lists its packages (NumPy, SciPy, HiGHS) without version pins, so record the
actual installed versions alongside any measurement you keep.

Notes for reading test output:

- `tests/known_failures.tsv`, `tests/qp_known_failures.tsv` and
  `tests/methods_findings.tsv` are tracked expected-failure and finding lists.
  They are part of the honest state of the solver, not hidden defects.
- Some benchmark scripts accept a second solver binary for side-by-side
  comparison; each script's docstring says exactly what it does and does not
  claim.
- CI timings are correctness signals, not benchmark figures.

## 5. Directory map

The maintained directory map is [docs/repository-map.md](docs/repository-map.md):
it lists the source, test, evidence and tooling directories, the build and
test entry points, and which
check to run for which kind of change. Read that first; it is updated in the
same commits that move files.

Orientation in four lines:

- `src/` is the C++17 solver core and the only code built into the evaluated
  engine; `tests/` and `benchmarks/` hold the shell, C++ and Python checks.
- `docs/` holds the written documentation; `results/` holds the committed
  evidence ledgers and is treated as read-only published records.
- `reproduction/` and `ci/` hold the pinned-corpus reproduction package and
  the per-commit CI gate.
- `scripts/dev.sh` is the single entry point for build, test, pinned-build and
  smoke tasks.

Retired material (the old single-file fallback engine and the old Python
`milp/` and `qp/` prototypes) lives on the
`archive/legacy-engines-2026-10-08` branch, not in main.

## 6. Where common edits live

| You want to change | Look at | Notes |
| --- | --- | --- |
| Model reading (MPS input) | `src/mps.cpp` | Python reference parsers under `parsers/`; semantics regressions in `benchmarks/mps_semantics_check.py` |
| Presolve | `src/presolve.cpp` | Enable with `--presolve`; replay logs via `--presolve-log`; design notes in `docs/verified_presolve.md`; tests in `benchmarks/presolve_tests.py`, `benchmarks/presolve_dominated_tests.py`, `benchmarks/presolve_replay_check.py` |
| Command-line parameters | `src/main.cpp` | The usage string and option parsing live together here; update both and extend `tests/cli_test.sh` |
| LP methods | `src/simplex.cpp`, `src/dual.cpp`, `src/ipm.cpp` | Certificate checks in `src/kkt_gate.cpp`, `src/nonoptimal_certificate.cpp` |
| Mixed-integer search | `src/milp.cpp` | Diagnostic flags `--audit-prop`, `--no-prop-prune`; suites under `tests/` (including `tests/milp_edge/` fixtures) and `benchmarks/milp_tests.py`; `tools/milp_check.py` is the random-MILP generator and HiGHS comparison |
| Benchmarks and checks | `benchmarks/`, `ci/` | New measured claims need the evidence workflow in section 7, not just a passing run |
| Public interface surface | `src/taral.hpp`, `docs/interface.md` | Keep the documented and actual behavior in sync |

`docs/foundation.md` carries the maintained module-responsibility table; treat
it as the authoritative map of the engine's internals.

## 7. How the gate and evidence workflow hangs together

The repository's rule is that every public claim traces to a retained
measurement. The moving parts:

1. **Measurement runners** (`benchmarks/`, `ci/gate.py`,
   `reproduction/reproduce.sh`) execute the solver against pinned corpora and
   write per-case ledgers. Corpora are hash-pinned
   (`benchmarks/netlib_manifest.sha256`,
   `reproduction/MANIFEST_corpus_sha256.txt`) so a rerun measures the same
   bytes.
2. **Independent checkers** (for example `benchmarks/orig_check.py`,
   `benchmarks/simplex_certificate_check.py`) re-verify reported answers
   against the original model files, separate from the solver's own
   arithmetic.
3. **Ledgers and provenance** land under `results/`, one directory per
   measurement, keeping the per-case table, run conditions and source
   snapshot together.
4. **Documents** (`docs/measured-results.md`, `docs/limitations.md`, the
   README requirement table) summarize ledgers into claims, with failures and
   unresolved cases retained alongside passes.
5. **The CI smoke** (`ci/smoke.sh`) re-checks a small pinned slice on every
   push so regressions surface before they reach a measured claim.

If you add or change solver behavior, the expected trail is: code, a test
next to the change, a ledger if the change is measurable, and wording that
matches what the ledger actually shows. Time-limited, approximate and
reference-missing outcomes are recorded as such, never silently dropped or
rounded up.

## 8. Making and checking your own edits

The license in `LICENSE` controls who may access this repository and for what
purpose; read it first. Within that scope, the intended evaluation loop is:

1. Edit on a local branch or copy. Keep one logical change per unit so its
   effect stays attributable.
2. Build warning-free with the strict flags from section 2.
3. Run `scripts/dev.sh test` (unit and CLI tests via ctest) and
   `scripts/dev.sh smoke` as the fast signal, then `bash ci/smoke.sh`
   before drawing conclusions.
4. Add or extend a test next to your change: shell cases in `tests/`,
   benchmark checks in `benchmarks/`.
5. If you measure anything, keep the run conditions with the numbers: source
   snapshot, compiler flags, machine, limits. A number without its conditions
   is not usable evidence.
6. Report findings or proposed changes through the repository's issue
   tracker, including the exact commit you evaluated.

## 9. Further reading

| Document | Topic |
| --- | --- |
| `docs/repository-map.md` | Maintained directory map, entry points, which check to run when |
| `README.md` | Requirement-by-requirement measured state with evidence links |
| `docs/interface.md` | Command-line and C++ interface reference and scope |
| `docs/measured-results.md` | Full measurement tables and verification rules |
| `docs/limitations.md` | Known limits and investigated discrepancies |
| `docs/foundation.md` | Engine module layout and dependency inspection |
| `docs/verified_presolve.md` | Presolve design and audit mechanism |
| `docs/deterministic_solve.md` | Deterministic work-budget behavior |
| `reproduction/README_REPRO.md` | Reproduction protocol, terminology and offline options |
