# TARAL-LP reproduction (C++ engine)

## Reading the technical terms

LP means linear programming; MILP means mixed-integer linear programming; QP means quadratic programming. IPM is an interior-point method. PDHG is primal-dual hybrid gradient, the approximate GPU method. KKT (Karush-Kuhn-Tucker) checks test feasibility and optimality conditions. A "gate" is the stated validation rule, not an exact-arithmetic proof. "Strict" means the separately stated tighter tolerance; a non-strict pass meets the ordinary rule but not that tighter check. A ledger is a per-case result table. Source hashes, file paths and command flags are retained only so engineers can reproduce a measurement. Historical measurements are not current-source claims.

JSON is the machine-readable result format. Warm starts reuse a previous solution or simplex basis; cold starts do not. fp64 is double-precision floating point; FMA is fused multiply-add.


From the root of a fresh clone, with Python dependencies installed and g++/git available:

```bash
python3 -m pip install -r reproduction/requirements.txt
bash reproduction/reproduce.sh
```

The script fetches the pinned Netlib corpus, verifies all file hashes, builds a fresh copy of the evaluated `src/` engine (`g++ -O3 -march=native -std=c++17`, the flags used for the committed ledgers), and runs the 93-case gate plus TRUSS at 60 seconds per case. Compile time is reported separately from solver time. It does not build the retired single-file fallback engine (kept on the `archive/legacy-engines-2026-10-08` branch). Historical result for the evaluated source snapshot recorded in that folder: 90/93 passes at 60 s (see `results/cpp_a9e8218_60s/`); this is not a current-source expected count; counts can vary by machine speed on slow cases.

Offline/subset check (the directory must contain the full corpus for manifest verification):

```bash
CORPUS=/absolute/path/to/mps_files bash reproduction/reproduce.sh afiro,e226
```

Docker uses the entire repository as its build context so canonical files are included:

```bash
docker build -f reproduction/Dockerfile -t taral .
docker run --rm taral
```

Docker has not been built or run in this environment. The reproduction folder is no longer a separately copyable package: the one-command flow requires this repository's `src/` and `benchmarks/` folders.

The engine uses only the C++ standard library. `benchmarks/netlib_gate.py` and `benchmarks/orig_check.py` are checking tooling (NumPy and HiGHS), not the engine solve path: HiGHS solves each original MPS file for the reference answer, and the independent checker re-reads the original file to verify every returned solution. Pass tolerances are fixed in `benchmarks/netlib_gate.py`. Output reports CPU only, no GPU. Requirements are unpinned; retain the tested environment versions with any submitted measurement.
