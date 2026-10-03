# CPU build and smoke checks

## Reading the technical terms

LP means linear programming; MILP means mixed-integer linear programming; QP means quadratic programming. IPM is an interior-point method. PDHG is primal-dual hybrid gradient, the approximate GPU method. KKT (Karush-Kuhn-Tucker) checks test feasibility and optimality conditions. A "gate" is the stated validation rule, not an exact-arithmetic proof. "Strict" means the separately stated tighter tolerance; a non-strict pass meets the ordinary rule but not that tighter check. A ledger is a per-case result table. Source hashes, file paths and command flags are retained only so engineers can reproduce a measurement. Historical measurements are not current-source claims.

JSON is the machine-readable result format. Warm starts reuse a previous solution or simplex basis; cold starts do not. fp64 is double-precision floating point; FMA is fused multiply-add.


From the repository root:

```sh
./scripts/build_and_smoke.sh
```

Requires Bash, CMake 3.16+, a C++17 compiler, and Python 3.10+ with `venv`
and `pip`. The script builds `build/taral` in Release mode, creates an isolated
`build/smoke-venv`, and installs the exact NumPy and HiGHS versions in
`requirements-smoke.txt`. The first dependency install needs internet access.
No system Python packages are changed. No CUDA, Netlib corpus, or SciPy is needed.

Overrides: `BUILD_DIR=/absolute/path`, `PYTHON=python3.11`, and `BUILD_JOBS=4`.
Optional runner arguments, such as `--seed 1`, follow the script name.

Build only, with no Python dependencies:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 2
```

Run again without installing anything:

```sh
build/smoke-venv/bin/python benchmarks/smoke.py --engine build/taral --out build/smoke
```

## Coverage and pass rule

The default seed is 57463. Existing benchmark generators and checks are reused:

- All 12 MPS semantic cases, including fixed/free format, objective constants,
  integer markers, and the quadratic objective convention.
- One existing LP fixture solved using simplex, dual simplex, and IPM. Returned
  column names, finite values, objective, rows, and bounds are independently
  checked by `orig_check.py`; HiGHS supplies the reference objective.
- Eight seeded MILP cases, one per existing suite, using one worker.
- Eight seeded QP cases, one per existing suite, plus the existing tiny
  quadratic-convention check, using one worker.
- The five existing Williams refinery scenarios, including the base published
  objective and agreement with HiGHS. This is not the illustrative refinery
  script, which imports a Python engine absent from this source snapshot.

A nonzero subprocess exit, timeout, missing benchmark summary, wrong result,
or unexpected undecided case fails the smoke command. The deliberately
node-limited `h_limits` MILP case may be undecided: the existing benchmark still
checks its incumbent and bound against HiGHS. This is recorded, not counted as
an optimal solve. Other undecided MILP/QP results fail the smoke run.

`build/smoke/summary.json` records UTC start time, seed, suite sizes, versions,
platform, binary/source/benchmark hashes, available source commit, available
CMake compiler/build settings, verdicts, and MILP/QP tallies. Per-suite stdout
and stderr are saved beside it. Source hashes identify actual solver contents
including uncommitted changes; a commit alone is not treated as provenance.
The build and logs are local artifacts, not intended to be committed.

This is a fast CPU regression check, not the full randomized or Netlib gate.
It does not establish performance, GPU correctness, or general numerical
reliability. No solver or existing benchmark behavior is changed.
