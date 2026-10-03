# QPLIB widening measurement

19 continuous convex-QP candidates. Engine requested cap 120s; HiGHS requested cap 300s. Input download and expanded-MPS SHA-256 values are recorded per row. Reported snapshot `58f77af` is independently bound by its source-file fingerprint. Build: g++13.3.0, `g++ -O2 -std=c++17 src/*.cpp -o taral`; reference HiGHS1.15.1; Kaggle Xeon2.20GHz,4 logical cores,31GB RAM. Caps are soft: engine ordinary timeout up to 140.68s and reference up to 363.74s. Objective gate is 1e-6*(1+|ref|); relative row/bound/integrality checks each 1e-6.

- Five returned objectives matched optimal HiGHS references under the reported comparison tolerance.
- Six feasible points were reported `optimal` by TARAL while HiGHS timed out or returned a solve error. No optimal reference was available for these six; feasibility and an `optimal` flag are not an independent optimality certificate. On 8845, 8495 and 8515, TARAL's objective is lower than HiGHS's time-limited incumbent.
- Four TARAL time limits, with HiGHS also time-limited. Two TARAL numerical failures, with HiGHS solve errors.
- 8991: the raw gate label `WRONG` was adjudicated using a local reproduction. Positive curvature (Gershgorin lower bound 2.0) and interior stationarity (gradient 1.12e-12) place the reproduced engine point numerically at the unique global minimum. The recorded HiGHS reference stopped short by 5.97e-6 absolute. This is not an established engine wrong answer; the raw CSV remains unchanged. See [diagnostic, inputs and limits](adjudication_8991.md).
- 9008 produced no JSON after 480.17s despite the requested 120s engine cap. Several ordinary `time_limit` rows also exceed 120s (10038:140.68s,8500:124.83s,8547:121.89s,8602:120.30s). A requested cap is not a guaranteed hard wall.

[Raw CSV](ledger.csv) · [Recorded provenance](provenance.json) · [Limitations](../../docs/limitations.md#qplib-measurements)
