# Repository map

Where things live, how to build, and how to run each kind of check. Paths are relative to the repository root.

## Build and run

| Task | Command |
| --- | --- |
| Build with warnings (CMake, Release) | `scripts/dev.sh build` -> `build/taral` |
| Build and run all unit and CLI tests | `scripts/dev.sh test` (ctest) |
| Exact measurement build used for the committed ledgers | `scripts/dev.sh pinned` -> `./taral` (`g++ -O3 -march=native -std=c++17`) |
| CPU-only build plus small deterministic checks | `scripts/dev.sh smoke` (see [benchmarks/SMOKE.md](../benchmarks/SMOKE.md)) |
| Pinned Netlib reproduction | [reproduction/README_REPRO.md](../reproduction/README_REPRO.md) |

CMake options (both default OFF, so a plain configure is unchanged): `-DTARAL_WARNINGS=ON` adds `-Wall -Wextra -Wpedantic`; `-DTARAL_BUILD_TESTS=ON` builds the C++ unit tests and registers them with `ctest`. The single executable target is `taral`. A plain `cmake` configure without `-DCMAKE_BUILD_TYPE` applies no optimisation flags; pass `-DCMAKE_BUILD_TYPE=Release` for speed.

## Directories

| Directory | What it holds |
| --- | --- |
| `src/` | The solver: MPS parser, simplex (primal and dual), interior point, MILP branch and bound, presolve, certificates, KKT gate, CLI. One CMake target, `taral`. [API header](../src/taral.hpp). |
| `tests/` | Shell and C++ tests for the CLI and gates (`cli_test.sh`, `parse_deadline_test.sh`, `presolve_work_limit_test.sh`, `kkt_gate_tests.cpp`, `ipm_gate/`), adversarial runners (`run_adversarial*.sh`), seeds (`seeds.txt`, `qp_seeds.txt`), known-failure tables, small fixtures and reproducer models, `milp_edge/` edge-case MILP models. |
| `tools/adv/`, `tools/advqp/` | Adversarial test harnesses and independent oracle checkers for LP/MILP and QP. The oracle MPS parser (`tools/adv/mpsio.py`) shares no code with `src/mps.cpp` on purpose. `tools/milp_check.py` is the random-MILP generator and HiGHS comparison. `tools/diag_qp_ipm/` holds QP diagnostics. |
| `benchmarks/` | Benchmark gates and ledgers (`netlib_gate.py`, `milp_ledger.py`, `qp_maros_meszaros.py`, ...), per-feature checkers (`*_check.py`) and tests (`*_tests.py`, `*_api_tests.cpp`), plotting, fixtures, and completed CPU-wave CSVs in `benchmarks/results/`. |
| `ci/` | CI gate: `smoke.sh` runs the pinned build, a ten-case Netlib gate for primal and dual, and the edge suites; `gate.py` and `mpsio.py` support it. |
| `reproduction/` | Pinned Netlib reproduction (Dockerfile, corpus hash manifest, `reproduce.sh`). |
| `results/` | Committed, manifest-pinned evidence ledgers. Read-only: treat as published records, not working files. |
| `docs/` | Method notes, measured results, limitations, interface and per-feature documents. Start from [measured-results.md](measured-results.md). |
| `examples/` | Synthetic refinery LP/MILP examples (no operational data). |
| `gpu/` | CUDA PDHG prototype and its measurements. Separate from the CMake build. |
| `experiments/` | Non-shipping prototypes with their own build scripts (`r2_miqp/`). |
| `engine/`, `parsers/` | Historical Python reference solver and MPS readers. Used by `examples/` and a few checks in `benchmarks/`; not part of the C++ engine. |
| `scripts/` | `dev.sh` (entry point above) and `build_and_smoke.sh`. |

Retired material (the single-file fallback engine and its snapshots, plus the old Python `milp/` and `qp/` prototypes) is kept on the `archive/legacy-engines-2026-10-08` branch.

## Which check to run

| You changed | Run |
| --- | --- |
| Anything in `src/` | `scripts/dev.sh test`, then `scripts/dev.sh smoke` |
| Parsing (`src/mps.cpp`) | `benchmarks/mps_semantics_check.py --engine <taral>` and `tests/cli_test.sh <taral>` |
| LP, MILP or QP behaviour | `tests/run_adversarial.sh` / `tests/run_adversarial_qp.sh` ([adversarial.md](adversarial.md), [adversarial_qp.md](adversarial_qp.md)) |
| Anything affecting a published number | the matching gate or ledger script in `benchmarks/`, compared against `results/` |

The adversarial runners and ledger scripts need Python 3 with NumPy, SciPy and HiGHS (`highspy`); see `benchmarks/requirements-smoke.txt` and `ci/requirements.txt`.
