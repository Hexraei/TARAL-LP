# TARAL-LP reproduction (C++ engine)
bash reproduce.sh            # clones the pinned Netlib corpus, verifies file hashes, cold-builds, runs the 93-case gate + truss at 60 s
CORPUS=/path/to/mps bash reproduce.sh afiro,e226    # offline / subset smoke test
docker build -t taral . && docker run --rm taral     # Dockerfile provided, NOT yet built or run (no Docker in my workspace)
Engine: src/taral.cpp, std only. tools/ and data/ are checker only (scipy and HiGHS reference). Backend printed: CPU only, no GPU.
Cold compile time is printed separately from per-case solver wall time. Pass rule is fixed in tools/gate_cpp.py; feasibility is recomputed from the returned vectors on the ORIGINAL MPS.
Versions are printed at run start; requirements.txt is unpinned, pin from the Kaggle record (Python 3.12.13, scipy 1.16.3) before submission.
