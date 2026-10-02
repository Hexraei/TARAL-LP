# TARAL-LP reproduction (C++ engine)

From the root of a fresh clone, with Python dependencies installed and g++/git available:

```bash
python3 -m pip install -r reproduction/requirements.txt
bash reproduction/reproduce.sh
```

The script fetches the pinned Netlib corpus, verifies all file hashes, cold-builds the canonical bundled engine in `cpp-engine/engine/taral.cpp`, and runs the 93-case gate plus TRUSS at 60 seconds per case. Compile time is reported separately from solver time. It does not build the separate modular `src/` engine or change that engine's README results.

Offline/subset check (the directory must contain the full corpus for manifest verification):

```bash
CORPUS=/absolute/path/to/mps_files bash reproduction/reproduce.sh afiro,e226
```

Docker uses the entire repository as its build context so canonical files are included:

```bash
docker build -f reproduction/Dockerfile -t taral .
docker run --rm taral
```

Docker has not been built or run in this environment. The reproduction folder is no longer a separately copyable package: the one-command flow requires this repository's `cpp-engine/` folder. Source, tools, reference JSON and fixtures have one canonical copy there.

The engine uses only the C++ standard library. `cpp-engine/tools/` and `cpp-engine/data/` are diagnostic/reference tooling (NumPy, SciPy and HiGHS), not the engine solve path. Pass tolerances remain fixed in `cpp-engine/tools/gate_cpp.py`, checking original-MPS feasibility/objective. Output reports CPU only, no GPU. Requirements are unpinned; retain the tested environment versions with any submitted measurement.
