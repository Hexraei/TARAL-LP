# TARAL-LP

A C++17 optimization solver core for refinery planning and related sparse industrial models. The current engine implements CPU LP, MILP and convex QP paths without linking an existing solver library. HiGHS is used for benchmark reference answers, not inside the solve path.

This is an experimental solver, not a production replacement for established industrial solvers. The strongest measured result is the CPU dual-simplex Netlib gate: 91/93 passes at 60 seconds per case, with no wrong answers found under that protocol. The default primal result is separate: 90/93.

## Problem statement coverage

The checklist separates implemented features from measured coverage and unfinished industrial requirements. DONE applies only to the named scope, not to the whole solver.

| Requirement | Status and evidence | What remains |
| --- | --- | --- |
| A sovereign solver core built from mathematical foundations, not an existing solver library | DONE for the current C++ engine: standard-library-only solve path in [`src/`](src/) | Broader industrial validation and performance work |
| A core with a basic API or CLI, rather than a modeling environment or GUI | DONE: C++ interfaces in [`src/taral.hpp`](src/taral.hpp), MPS input and CLI below | Stable versioned API and broader integration testing |
| LP as an initial focus | DONE for measured Netlib coverage: 91/93 dual and 90/93 primal under the stated gate | Unresolved cases and larger industrial models |
| MILP as an initial focus | IN PROGRESS: branch-and-bound implemented; selected MIPLIB measurements below | Faster search and broader difficult-MILP coverage |
| QP as an initial focus | IN PROGRESS: convex-QP interior-point path and [Maros-Meszaros ledger](results/qp_maros_meszaros/ledger.csv) | Reference gaps, numerical failures and uncertified cases; nonconvex QP is rejected |
| Modular extension to MIQP, NLP and MINLP | IN PROGRESS: solver modules are separate in [`src/`](src/) | These problem classes are not implemented or measured |
| Revised simplex for continuous problems | DONE: primal and dual CPU simplex, measured on Netlib | Scaling and speed improvements |
| Interior-point methods for continuous problems | IN PROGRESS: CPU interior-point path in [`src/`](src/), [LP results](results/ipm_netlib/) and QP results | Broader convergence and performance validation |
| Branch-and-bound for MILP | DONE as an implementation: [`src/milp.cpp`](src/milp.cpp), randomized tests and MIPLIB measurements | Industrial-scale search performance |
| Branch-and-cut, cutting planes, presolve, heuristics and node selection | IN PROGRESS: baseline MILP search exists | Full requirement coverage and separate measured evidence for each search technique |
| Sparse matrices and efficient numerical linear algebra | DONE as an implementation: sparse columns and sparse LU in [`src/`](src/) | Further memory/scaling measurements and factorization tuning |
| Multi-core parallelization | IN PROGRESS | Current headline engine runs single-threaded; parallel solve performance is not established |
| GPU acceleration where it brings a measured benefit | IN PROGRESS: approximate CUDA PDHG measurements below | Equal-accuracy crossover evidence and stricter original-model feasibility |
| Numerical stability and reliable convergence | IN PROGRESS: original-model checks, strict residual results and known unresolved cases below | Full robustness across ill-conditioned and difficult industrial models |
| Thousands to millions of variables, sparse and highly constrained industrial models | IN PROGRESS | Selected benchmarks do not establish consistent industrial-scale performance |
| Degenerate models, ill-conditioned matrices and weak MILP relaxations | IN PROGRESS: Netlib cases and [`benchmarks/milp_tests.py`](benchmarks/milp_tests.py) exercise these areas | Broader independently checked stress coverage and practical solve times |
| Refinery scheduling, crude blending and process optimization | IN PROGRESS: synthetic refinery examples in [`examples/`](examples/) | Operational data, scheduling coverage and measured industrial outcomes |
| Production planning, logistics, power dispatch, transportation and supply chain | IN PROGRESS | Validated case studies for each application, not only generic benchmark models |
| Netlib, MIPLIB and Mittelmann benchmarks compared with an established solver | IN PROGRESS: Netlib and selected MIPLIB/Kennington comparisons with HiGHS below | Completed Mittelmann evidence and broader benchmark coverage |
| QPLIB where applicable | IN PROGRESS: 2/2 selected CPU guard measurements below | Convex-QP solve coverage on applicable QPLIB instances |
| Public benchmarks plus open-literature refinery/blending/planning/supply-chain case studies | IN PROGRESS: public benchmark ledgers and synthetic examples | Traceable open-literature cases and all named application areas |
| A transparent, extensible foundation | DONE for inspectable source and retained result ledgers | Stable proof interfaces, broader reproducibility and industrial validation |

## What sets the project apart

### What we are building toward: a solver that proves its answers

The goal is that a returned LP answer comes with evidence a separate checker can verify: KKT evidence for an optimum, a Farkas proof for infeasibility, or a feasible anchor and improving ray for unboundedness. A checker should not need to trust the solver's own success flag.

Current main uses independent original-model benchmark checks. The `taral-ai/simplex-certificates` and `taral-ai/infeasible-unbounded-certificates` branches contain certificate development; they are unmerged. Certificate export is development work, not a promise that every solve on main already ships a proof. Limits and numerical failures remain unresolved outcomes. MILP global proofs and nonconvex-QP certificates are outside this claim.

### What we are building toward: a solver built for refinery problems

The model structure supports sparse balances, row ranges, bounds and mixed-integer decisions. Current checked-in refinery examples are synthetic, not field data. The unmerged `taral-ai/refinery-stress` branch has twelve refinery-shaped synthetic fixtures checked against HiGHS; those tests do not establish industrial readiness and are not yet part of main.

The target is a solver core that fits refinery models, not a GUI wrapped around another solver. Measured refinery savings, live plant integration and production scheduling remain unverified.

## Measured CPU results

### Netlib LP

An independent Kaggle run checked pinned `src/` at `442ca16` against live HiGHS references on the original MPS files. Local ledgers retain the primal and dual protocols separately.

| Protocol | Outcome | Evidence |
| --- | --- | --- |
| Default primal, 60 seconds per case | 90/93 passes; 0 wrong answers under this gate | [`results/cpp_a9e8218_60s/`](results/cpp_a9e8218_60s/) |
| Stricter primal check | 89/90 passing cases meet row/bound violation and relative objective error <= 1e-8 | Same primal ledger and independent Kaggle check |
| Dual simplex, 60 seconds per case | 91/93 passes; 0 wrong answers; strict check 90/91 | [`results/cpp_0bad060_dual_60s/`](results/cpp_0bad060_dual_60s/) |
| Separate extended primal protocol, 300 seconds per case | 92/93 passes; 0 wrong answers | [`results/cpp_a9e8218_ext300/`](results/cpp_a9e8218_ext300/) |

The 93-case denominator excludes PILOT.WE and PILOT4 under the fixed local 60-second reference rule. Their separate Kaggle/extended results do not change that denominator. TRUSS passes as an extra case. DFL001 reaches the default-primal time limit; dual finishes in about 36 seconds on Kaggle. GREENBEA clears the fixed gate but misses the stricter check, with measured row violation 1.46e-8 on a cancellation-heavy row.

### CPU benchmark wave - October 2, 2026

All listed runs completed in a Kaggle CPU notebook using pinned engine source `442ca16`. "Measured" counts completed runs, not optimal solves. TARAL had a 30-second cap per case.

| Benchmark selection | Completed measurements | Outcome | Ledger |
| --- | --- | --- | --- |
| MIPLIB | 12/12 measured | Six optimal objectives matched HiGHS; six reached the 30-second time limit | [`miplib_ledger.csv`](benchmarks/results/miplib_ledger.csv) |
| Kennington LP | 12/12 measured | Six optimal objectives matched HiGHS; six reached the 30-second time limit | [`kennington_ledger.csv`](benchmarks/results/kennington_ledger.csv) |
| QPLIB | 2/2 measured | Both were rejected by the nonconvex-QP guard, which reported that Q was not positive semidefinite in the minimization sense | [`qplib_ledger.csv`](benchmarks/results/qplib_ledger.csv) |

MIPLIB objective matches: p0033, p0201, mod008, stein27, misc03 and lseu. Kennington objective matches: ken-07, cre-a, cre-c, osa-07, osa-14 and osa-30. The QPLIB guard measurements are QPLIB_0018 and QPLIB_0343; these are rejection checks, not convex-QP solves. For stein45, HiGHS also reached its separate 60-second cap, leaving no optimal reference.

These are selections, not full-library coverage. The wave CSVs record objective/status comparisons, not a separate point-feasibility or certificate check. No speed advantage is claimed from them.

A local check of 3,600 small LPs with known optimal/infeasible/unbounded status found zero false verdicts.

### Other retained CPU measurements

- Small MIPLIB selection: 26/26 runs recorded; 19 solved at the 300-second cap, seven reached that cap. The separate 30-second protocol solved fourteen. No wrong answer was found under those checks. See [`results/milp_miplib/`](results/milp_miplib/).
- Convex-QP Maros-Meszaros set: 99/99 measured; 71 passed the stated rule. Remaining rows include five objective discrepancies, eighteen without a usable HiGHS reference, four numerical failures and one nonconvex rejection. Four objective discrepancies were traced to reference accuracy/parsing; the fifth remains uncertified. See the [unchanged ledger and notes](results/qp_maros_meszaros/NOTES.md).
- Local small-MILP holdout at `442ca16`: 300/300 measured with seeds 26120, 31 and 47. 296 cleared the mismatch flag, including four separately recorded reference-objective discrepancies with feasible integral TARAL points. Three original inputs had undeclared-column parser errors; corrected inputs matched HiGHS. One flag was a reference-side presolve discrepancy. This was local, not Kaggle-gated or industrial-scale validation.

## Approximate GPU LP results

CUDA PDHG is separate from the exact-simplex headline. The Kaggle T4 check used `gpu/pdhg.cu` from `442ca16`, `sm_75`, tolerance `1e-6`, three repeats per case, live SciPy/HiGHS references and original-model point checks.

| Case | GPU solve time | Same-run CPU context |
| --- | --- | --- |
| DFL001 | 0.60-0.71 s | Primal reached 90-second cap; dual 31.8 s |
| PILOT87 | 26.8-26.9 s | CPU PDHG reached 120-second cap unconverged; dual 27.6-28.2 s |
| FIT2P | 6.3-6.5 s | Dual about 5.0 s |
| Synthetic transport_316x316, seed 1 | 3.7-3.8 s | Dual about 0.5 s |

All twelve GPU runs were measured. None met the strict `1e-6` original-model row/bound/objective gate: relative row violations ranged from `1.2e-5` to `3.4e-3`. Relative objective errors ranged from about `2e-11` to `1.4e-6`. Setup adds 0.35-0.77 seconds separately. Self-reported dual residuals and gaps are not independent certificates.

These selected cases show fast approximate answers, not GPU wins at equal accuracy. FIT2P and transport favor CPU dual simplex. See [`results/gpu_pdhg/`](results/gpu_pdhg/).

A separate Kaggle T4 synthetic scaling table shows why acceleration is model-dependent. CPU dual is faster on every listed transport size through 1,000,000 variables (8.34 seconds at that size, versus about 48.15 seconds GPU solve time plus setup). For packing at 10,000 variables, CPU dual reaches its 90-second cap while GPU returns an approximate answer in 0.20-0.35 seconds of solve time. Larger listed packing cases show the same runtime pattern, but no packing row meets the table's strict `1e-6` gate, and objective error against HiGHS is unavailable from 100,000 variables onward. This is a speed observation at different achieved accuracy, not proof of faster correct solves. See [`gpu_scale_T4_table.csv`](benchmarks/results/gpu_scale_T4_table.csv).

## Build and use

```bash
g++ -O3 -march=native -std=c++17 -o taral src/*.cpp
./taral path/to/model.mps --time-limit 60 --sol out.sol --json out.json
./taral path/to/model.mps --method dual --time-limit 60 --json out.json
```

The engine has no external solver or linear-algebra dependency. The C++ API is in [`src/taral.hpp`](src/taral.hpp). JSON reports status, objective, iterations, wall time and message. A solution file reports original variable values for an optimal LP. Time limits, numerical failures, unsupported inputs and parse errors must not be read as solved cases.

## Verification and reproduction

The fixed Netlib pass rule requires optimal status from both engines, objective error <= `max(1e-6, 1e-7 * abs(HiGHS objective))`, independently recomputed objective agreement, relative row violation <= `1e-6` and bound violation <= `1e-6`. Strict results use the tighter check separately. Each other benchmark retains its own ledger/protocol; “measured” does not replace its solve-quality checks.

Follow [`reproduction/README_REPRO.md`](reproduction/README_REPRO.md) for the C++ Netlib build, pinned corpus and 93-case gate. Local headline hardware: Intel Core 7 240H, 10 cores/16 threads, 15 GiB RAM, Ubuntu 24.04, Linux 6.17, g++ 13.3, `-O3 -march=native -std=c++17`. Cases ran one at a time on a single-threaded engine. Local and Kaggle timing records are separate; near-cap outcomes depend on host speed.

The historical six-variable refinery example retains a hypothetical objective of about $1.19 million/day in assumed units. It is not a measured margin or savings result. Run it with `python -m examples.illustrative_refinery` after installing NumPy and SciPy.

## Repository guide

- [`src/`](src/): current C++ solver, parser and CLI.
- [`benchmarks/`](benchmarks/): reference checks and benchmark drivers; [`benchmarks/results/`](benchmarks/results/) holds the completed CPU wave CSVs.
- [`results/`](results/): retained C++ and historical ledgers. Historical files are not current-engine claims.
- [`gpu/`](gpu/): CUDA kernels and approximate PDHG prototype.
- [`reproduction/`](reproduction/): pinned C++ Netlib reproduction.
- [`examples/`](examples/): synthetic refinery examples, without operational-data or savings validation.
- [`cpp-engine/`](cpp-engine/): separate bundled reference engine; its results do not apply to `src/`.
- [`engine/`](engine/), [`parsers/`](parsers/), [`milp/`](milp/), [`qp/`](qp/): historical Python prototypes, separate from the current C++ engine.

## License and access

Proprietary. Only official competition judges and organizers may read this submission for evaluation. No copying, modification, redistribution, commercial use or AI/ML use; see [LICENSE](LICENSE) for full terms. The repository remains private until its owner chooses otherwise.
