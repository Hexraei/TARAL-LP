# TARAL-LP reproduction (C++ engine)

From the root of a fresh clone, with Python dependencies installed and g++/git available:

```bash
python3 -m pip install -r reproduction/requirements.txt
bash reproduction/reproduce.sh
```

The script fetches the pinned Netlib corpus, verifies all file hashes, cold-builds the headline `src/` engine (`g++ -O3 -march=native -std=c++17`, the flags used for the committed ledgers), and runs the 93-case gate plus TRUSS at 60 seconds per case. Compile time is reported separately from solver time. It does not build the separate single-file `cpp-engine/` fallback. Expected result on the committed engine: 90/93 passes at 60 s (see `results/cpp_a9e8218_60s/`); the exact count can vary by machine speed on slow cases.

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
