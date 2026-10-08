# Detailed measurements and retained context

## Reading the technical terms

LP means linear programming; MILP means mixed-integer linear programming; QP means quadratic programming. IPM is an interior-point method. PDHG is primal-dual hybrid gradient, the approximate GPU method. KKT (Karush-Kuhn-Tucker) checks test feasibility and optimality conditions. A "gate" is the stated validation rule, not an exact-arithmetic proof. "Strict" means the separately stated tighter tolerance; a non-strict pass meets the ordinary rule but not that tighter check. A ledger is a per-case result table. Source hashes, file paths and command flags are retained only so engineers can reproduce a measurement. Historical measurements are not current-source claims.

JSON is the machine-readable result format. Warm starts reuse a previous solution or simplex basis; cold starts do not. fp64 is double-precision floating point; FMA is fused multiply-add.


This is the long-form README retained during the October 3 layout revision. Historical tables and implementation notes refer to their stated snapshots, not an updated claim that every historical implementation note describes the current solver. [Current landing page](../README.md) · [Current limitations](limitations.md)

# TARAL-LP

A C++17 optimization solver core for refinery planning and related sparse industrial models. The current engine implements CPU LP, MILP and convex QP paths without linking an existing solver library. HiGHS is used for benchmark reference answers, not inside the solve path.

This is an experimental solver, not a production replacement for established industrial solvers. The post-merge CPU Netlib gate measured 93/93 passes in both primal and dual modes, with no wrong answers and 92/93 meeting the stricter check ([ledgers](../results/netlib_postmerge_58f77af/)).

Protocol: 60s per engine solve, PILOT.WE and PILOT4 included against HiGHS references capped at 300s; default `reproduction/reproduce.sh` excludes `pilot.we`/`pilot4` from passes as `reference_excluded_60s` (ceiling 91/93).

Project website: [taral-lp.vercel.app](https://taral-lp.vercel.app/) - interactive performance reports, per-case ledgers and measured limits. The full site is prepared for publication; the live URL currently shows the holding page.

## Problem statement coverage

The checklist separates implemented features from measured coverage and unfinished industrial requirements. DONE applies only to the named scope, not to the whole solver.

| Requirement | Status and evidence | What remains |
| --- | --- | --- |
| A sovereign solver core built from mathematical foundations, not an existing solver library | DONE for the current C++ engine: standard-library-only solve path in [`src/`](../src/) | Broader industrial validation and performance work |
| A core with a basic API or CLI, rather than a modeling environment or GUI | DONE: C++ interfaces in [`src/taral.hpp`](../src/taral.hpp), MPS input and CLI below | Stable versioned API and broader integration testing |
| LP as an initial focus | DONE for measured Netlib coverage: 93/93 in both modes under the validation including both PILOT.WE and PILOT4 ([protocol above](#taral-lp)) | Unresolved cases and larger industrial models |
| MILP as an initial focus | IN PROGRESS: branch-and-bound implemented; selected MIPLIB measurements below | Faster search and broader difficult-MILP coverage |
| QP as an initial focus | IN PROGRESS: convex-QP interior-point path and [Maros-Meszaros ledger](../results/qp_maros_meszaros/ledger.csv) | Reference gaps, numerical failures and uncertified cases; nonconvex QP is rejected |
| Modular extension to MIQP, NLP and MINLP | IN PROGRESS: solver modules are separate in [`src/`](../src/) | These problem classes are not implemented or measured |
| Revised simplex for continuous problems | DONE: primal and dual CPU simplex, measured on Netlib | Scaling and speed improvements |
| Interior-point methods for continuous problems | IN PROGRESS: CPU interior-point path in [`src/`](../src/), [LP results](../results/ipm_netlib/) and QP results | Broader convergence and performance validation |
| Branch-and-bound for MILP | DONE as an implementation: [`src/milp.cpp`](../src/milp.cpp), randomized tests and MIPLIB measurements | Industrial-scale search performance |
| Branch-and-cut, cutting planes, presolve, heuristics and node selection | IN PROGRESS: baseline MILP search exists | Full requirement coverage and separate measured evidence for each search technique |
| Sparse matrices and efficient numerical linear algebra | DONE as an implementation: sparse columns and sparse LU in [`src/`](../src/) | Further memory/scaling measurements and factorization tuning |
| Multi-core parallelization | IN PROGRESS | Current headline engine runs single-threaded; parallel solve performance is not established |
| GPU acceleration where it brings a measured benefit | IN PROGRESS: approximate CUDA PDHG and measured crossover below | No exact-answer crossover win; basis identification, dense-nucleus factorization and broader original-model checks remain |
| Numerical stability and reliable convergence | IN PROGRESS: original-model checks, strict residual results and known unresolved cases below | Full robustness across ill-conditioned and difficult industrial models |
| Thousands to millions of variables, sparse and highly constrained industrial models | IN PROGRESS | Selected benchmarks do not establish consistent industrial-scale performance |
| Degenerate models, ill-conditioned matrices and weak MILP relaxations | IN PROGRESS: Netlib cases and [`benchmarks/milp_tests.py`](../benchmarks/milp_tests.py) exercise these areas | Broader independently checked stress coverage and practical solve times |
| Refinery scheduling, crude blending and process optimization | IN PROGRESS: synthetic refinery examples in [`examples/`](../examples/) | Operational data, scheduling coverage and measured industrial outcomes |
| Production planning, logistics, power dispatch, transportation and supply chain | IN PROGRESS | Validated case studies for each application, not only generic benchmark models |
| Netlib, MIPLIB and Mittelmann benchmarks compared with an established solver | IN PROGRESS: Netlib and selected MIPLIB/Kennington comparisons with HiGHS below | Completed Mittelmann evidence and broader benchmark coverage |
| QPLIB where applicable | IN PROGRESS: 2/2 selected CPU guard measurements below | Convex-QP solve coverage on applicable QPLIB instances |
| Public benchmarks plus open-literature refinery/blending/planning/supply-chain case studies | IN PROGRESS: public benchmark ledgers and synthetic examples | Traceable open-literature cases and all named application areas |
| A transparent, extensible foundation | DONE for inspectable source and retained result ledgers | Stable proof interfaces, broader reproducibility and industrial validation |

## What sets the project apart

### What we are building toward: a solver that proves its answers

The goal is that a returned LP answer comes with evidence a separate checker can verify: KKT evidence for an optimum, a Farkas proof for infeasibility, or a feasible anchor and improving ray for unboundedness. A checker should not need to trust the solver's own success flag.

The current solver uses independent original-model benchmark checks. Historical development snapshots contain optimality and infeasibility/unboundedness certificate work; those snapshot notes do not establish the current integration state. Certificate export is development work, not a promise that every solve already exports a verified proof. Limits and numerical failures remain unresolved outcomes. MILP global proofs and nonconvex-QP certificates are outside this claim.

### What we are building toward: a solver built for refinery problems

The model structure supports sparse balances, row ranges, bounds and mixed-integer decisions. Current checked-in refinery examples are synthetic, not field data. A historical development snapshot has twelve refinery-shaped synthetic fixtures checked against HiGHS; those tests do not establish industrial readiness and are not established here as part of the evaluated solver.

The target is a solver core that fits refinery models, not a GUI wrapped around another solver. Measured refinery savings, live plant integration and production scheduling remain unverified.

## Measured CPU results

### Netlib LP

The post-merge Kaggle gate on the evaluated source snapshot `58f77af` measured 93/93 in both modes, strict 92/93, with only GREENBEA outside the tighter tolerance ([raw ledgers and logs](../results/netlib_postmerge_58f77af/)); the counting protocol is stated above. Its `src/` is byte-identical to the evaluated source snapshot `6c621e0`.

The historical table below is unchanged: an independent Kaggle run checked pinned `src/` at `442ca16` against live HiGHS references on the original MPS files. Local ledgers retain the primal and dual protocols separately.

| Protocol | Outcome | Evidence |
| --- | --- | --- |
| Default primal, 60 seconds per case | 90/93 passes; 0 wrong answers under this gate | [`results/cpp_a9e8218_60s/`](../results/cpp_a9e8218_60s/) |
| Stricter primal check | 89/90 passing cases meet row/bound violation and relative objective error <= 1e-8 | Same primal ledger and independent Kaggle check |
| Dual simplex, 60 seconds per case | 91/93 passes; 0 wrong answers; strict check 90/91 | [`results/cpp_0bad060_dual_60s/`](../results/cpp_0bad060_dual_60s/) |
| Separate extended primal protocol, 300 seconds per case | 92/93 passes; 0 wrong answers | [`results/cpp_a9e8218_ext300/`](../results/cpp_a9e8218_ext300/) |

The local ledgers contain 94 rows: 93 denominator rows plus TRUSS as an extra case. PILOT.WE and PILOT4 are marked `reference_excluded_60s`, but their `in_denominator` fields remain true and the summaries retain a denominator of 93. They count as nonpasses in the published 60-second totals; do not drop them and change the denominator. The separate extended protocol does not alter the 60-second headline. DFL001 reaches the default-primal time limit; dual finishes in about 36 seconds on Kaggle. GREENBEA clears the fixed gate but misses the stricter check, with measured row violation 1.46e-8 on a cancellation-heavy row.

### CPU benchmark selection - October 2, 2026

All listed runs completed in a Kaggle CPU notebook using pinned engine source `442ca16`. "Measured" counts completed runs, not optimal solves. TARAL had a 30-second cap per case.

Measurement environment: Kaggle CPU notebook, no accelerator; Intel Xeon @ 2.20GHz, four logical CPUs (one socket, two cores, two threads per core), 30 GiB RAM cap; x86_64, Linux 6.18.48+, Ubuntu 22.04.5 LTS. Build: `g++ -O3 -march=native -std=c++17 -Wall -Wextra -Wpedantic -o taral src/*.cpp`.

| Benchmark selection | Completed measurements | Outcome | Ledger |
| --- | --- | --- | --- |
| MIPLIB | 12/12 measured | Six optimal objectives matched HiGHS; six reached the 30-second time limit | [`miplib_ledger.csv`](../benchmarks/results/miplib_ledger.csv) |
| Kennington LP | 12/12 measured | Six optimal objectives matched HiGHS; six reached the 30-second time limit | [`kennington_ledger.csv`](../benchmarks/results/kennington_ledger.csv) |
| QPLIB | 2/2 measured | Both were rejected by the nonconvex-QP guard, which reported that Q was not positive semidefinite in the minimization sense | [`qplib_ledger.csv`](../benchmarks/results/qplib_ledger.csv) |

MIPLIB objective matches: p0033, p0201, mod008, stein27, misc03 and lseu. Kennington objective matches: ken-07, cre-a, cre-c, osa-07, osa-14 and osa-30. The QPLIB guard measurements are QPLIB_0018 and QPLIB_0343; these are rejection checks, not convex-QP solves. For stein45, HiGHS also reached its separate 60-second cap, leaving no optimal reference.

These are selections, not full-library coverage. The wave CSVs record objective/status comparisons, not a separate point-feasibility or certificate check. No speed advantage is claimed from them.

A local check of 3,600 small LPs with known optimal/infeasible/unbounded status found zero false verdicts.

### Other retained CPU measurements

- Small MIPLIB selection: 26/26 runs recorded; 19 solved at the 300-second cap, seven reached that cap. The separate 30-second protocol solved fourteen. No wrong answer was found under those checks. See [`results/milp_miplib/`](../results/milp_miplib/).
- Convex-QP Maros-Meszaros set: 99/99 measured; 71 passed the stated rule. Remaining rows include five objective discrepancies, eighteen without a usable HiGHS reference, four numerical failures and one nonconvex rejection. Four objective discrepancies were traced to reference accuracy/parsing; the fifth remains uncertified. See the [unchanged ledger and notes](../results/qp_maros_meszaros/NOTES.md).
- Local small-MILP holdout at `442ca16`: 300/300 measured with seeds 26000+120, 31 and 47. 296 cleared the mismatch flag, including four separately recorded reference-objective discrepancies with feasible integral TARAL points. Three original inputs had undeclared-column parser errors; corrected inputs matched HiGHS. One flag was a reference-side presolve discrepancy. This was local, not Kaggle-gated or industrial-scale validation.

## Approximate GPU LP results

CUDA PDHG is separate from the exact-simplex headline. The Kaggle T4 check used `gpu/pdhg.cu` from the evaluated source snapshot `442ca16`, `sm_75`, tolerance `1e-6`, three repeats per case, live SciPy/HiGHS references and original-model point checks.

| Case | GPU solve time | Same-run CPU context |
| --- | --- | --- |
| DFL001 | 0.60-0.71 s | Primal reached 90-second cap; dual 31.8 s |
| PILOT87 | 26.8-26.9 s | CPU PDHG reached 120-second cap unconverged; dual 27.6-28.2 s |
| FIT2P | 6.3-6.5 s | Dual about 5.0 s |
| Synthetic transport_316x316, seed 1 | 3.7-3.8 s | Dual about 0.5 s |

All twelve GPU runs were measured. None met the strict `1e-6` original-model row/bound/objective gate: relative row violations ranged from `1.2e-5` to `3.4e-3`. Relative objective errors ranged from about `2e-11` to `1.4e-6`. Setup adds 0.35-0.77 seconds separately. Self-reported dual residuals and gaps are not independent certificates.

These selected cases show fast approximate answers, not GPU wins at equal accuracy. FIT2P and transport favor CPU dual simplex. See [`results/gpu_pdhg/`](../results/gpu_pdhg/).

A separate Kaggle T4 synthetic scaling table shows why acceleration is model-dependent. CPU dual is faster on every listed transport size through 1,000,000 variables (8.34 seconds at that size, versus about 48.15 seconds GPU solve time plus setup). For packing at 10,000 variables, CPU dual reaches its 90-second cap while GPU returns an approximate answer in 0.20-0.35 seconds of solve time. Larger listed packing cases show the same runtime pattern, but no packing row meets the table's strict `1e-6` gate, and objective error against HiGHS is unavailable from 100,000 variables onward. This is a speed observation at different achieved accuracy, not proof of faster correct solves. See [`gpu_scale_T4_table.csv`](../benchmarks/results/gpu_scale_T4_table.csv).

### Local GPU scaling and exact controls - October 2, 2026

These newer local runs are separate from the Kaggle T4 measurements above. Hardware: RTX 5050 Laptop GPU (8 GB, compute capability 12.0), Intel Core 7 240H and CUDA 12.9. The round-2 synthetic tables use source `5e93146`, a 120-second limit, three repeats per PDHG configuration and a single run per exact-simplex configuration. PDHG times below include setup and transfer. Finished PDHG vectors were rechecked on the original MPS files; no recheck failure was recorded. That is not an exact-answer certificate.

| Family and case | PDHG tolerance | GPU total | CPU PDHG, 10 threads | Exact dual, one thread | Evidence |
| --- | --- | --- | --- | --- | --- |
| Transport, 1,000,000 columns / 2,000,000 nonzeros | 1e-4 | 15.5 s | 74.6 s | 2.88 s | [transport 1e-4](../results/gpu_pdhg/synth_transport_tol1e-4.csv) |
| Transport, 5,000,000 columns / 10,000,000 nonzeros | 1e-6 | Time limit, about 120 s | Time limit, about 121 s | 19.4 s | [transport 1e-6](../results/gpu_pdhg/synth_transport_tol1e-6.csv) |
| Packing, 100,000 columns | 1e-4 | 0.743 s | 0.593 s | Time limit, 120 s | [packing 1e-4](../results/gpu_pdhg/synth_packing_tol1e-4.csv) |
| Packing, 1,000,000 columns | 1e-4 | 8.99 s | 9.31 s | Time limit, 120 s | [packing 1e-4](../results/gpu_pdhg/synth_packing_tol1e-4.csv) |
| Packing, 1,000,000 columns | 1e-6 | 30.3 s | 52.1 s | Time limit, 120 s | [packing 1e-6](../results/gpu_pdhg/synth_packing_tol1e-6.csv) |

Exact dual simplex beats every PDHG configuration on transport at every tested size and both tolerances. GPU PDHG beats CPU PDHG on some larger transport models, but that comparison misses the faster exact method. The largest transport model does not converge on GPU at 1e-6 within 120 seconds.

On packing, neither exact method finishes from 100,000 columns upward within the tested limit. PDHG is the only tested route returning a near-optimal answer there, not a certified exact answer. The 8.99 versus 9.31 second margin at one million columns and 1e-4 is small enough to sit within laptop run-to-run variation. Packing at one and two million columns has no usable reference objective, so objective accuracy is unknown there. See the [protocol and limits](../results/gpu_pdhg/NOTES.md).

The local GPU Netlib sweeps cover 94 files including TRUSS, not the CPU 93-case gate. At 1e-4, 92 report `near_optimal` and two reach the limit; at 1e-6, 86 report `near_optimal` and eight reach the limit. The tighter count changed near the cap between runs. At 1e-4, nine near-optimal rows have objective error above 1e-3 and the worst is about 4.8e-2 (FORPLAN). A scaled stopping tolerance is not a bound on objective error. Sources: [1e-4 sweep](../results/gpu_pdhg/netlib_94files_tol1e-4.csv), [1e-6 sweep](../results/gpu_pdhg/netlib_94files_tol1e-6.csv).

### GPU-to-simplex crossover: it does not pay yet

Crossover converts a PDHG point and row multipliers into a corner basis, then asks exact simplex to finish. Three alternating runs of each route on the same local laptop give DFL001 a median total of 20.4 seconds after handoff, versus 11.88 seconds for dual simplex without a previous solution. The previous-basis method needs 29,869 simplex iterations against 21,256 cold. These times include the PDHG process, not just its 0.59-second solve.

The point sits inside the optimal face rather than at a vertex. Only 3,507 of 6,071 basis slots are determined by variables off their bounds; the resulting basis has 2,585 dual infeasibilities and needs a 12,083-iteration dual phase 1. Tighter PDHG tolerances did not fix the basis. Packing with 100,000 columns has a nearly dense ~4,800 by ~4,800 factorization nucleus; neither the warm nor fresh-start method returns an exact answer in 120 seconds.

The factorization deadline now holds. The packing handoff that previously ran beyond 11.5 minutes returns `time_limit` at 60.07 seconds with a 60-second limit, instead of misreporting a singular basis or hanging. A warm basis gets a quarter of the limit to factor and falls back to the slack basis if it cannot. See [CROSSOVER.md, round 2](../results/gpu_pdhg/CROSSOVER.md) for the measurements and diagnosis.

Validation of the retained evaluated source snapshot `5e93146` gates still report dual 91/93 (strict 90/91) and primal 90/93 (strict 89/90), with zero wrong answers under that protocol: [dual run](../results/cpp_5e93146_dual_60s/summary.json), [primal run](../results/cpp_5e93146_primal_60s/summary.json). Faster approximate points have not become faster exact answers. The next crossover work is basis identification and factorization that can handle the dense nucleus.

## Build and use

```bash
g++ -O3 -march=native -std=c++17 -o taral src/*.cpp
./taral path/to/model.mps --time-limit 60 --sol out.sol --json out.json
./taral path/to/model.mps --method dual --time-limit 60 --json out.json
```

The engine has no external solver or linear-algebra dependency. The C++ API is in [`src/taral.hpp`](../src/taral.hpp). JSON reports status, objective, iterations, wall time and message. A solution file reports original variable values for an optimal LP. Time limits, numerical failures, unsupported inputs and parse errors must not be read as solved cases.

Exit codes: 0 definitive answer (optimal, infeasible, unbounded), 2 usage error, 3 model parse error, 4 time/iteration/node limit, 5 other failure (numerical, unsupported, nonconvex).

## Verification and reproduction

The fixed Netlib pass rule requires optimal status from both engines, objective error <= `max(1e-6, 1e-7 * abs(HiGHS objective))`, independently recomputed objective agreement, relative row violation <= `1e-6` and bound violation <= `1e-6`. Strict results use the tighter check separately. Each other benchmark retains its own ledger/protocol; “measured” does not replace its solve-quality checks.

Follow [`reproduction/README_REPRO.md`](../reproduction/README_REPRO.md) for the C++ Netlib build, pinned corpus and 93-case gate. Local headline hardware: Intel Core 7 240H, 10 cores/16 threads, 15 GiB RAM, Ubuntu 24.04, Linux 6.17, g++ 13.3, `-O3 -march=native -std=c++17`. Cases ran one at a time on a single-threaded engine. Local and Kaggle timing records are separate; near-cap outcomes depend on host speed.

The historical six-variable refinery example retains a hypothetical objective of about $1.19 million/day in assumed units. It is not a measured margin or savings result. Run it with `python -m examples.illustrative_refinery` after installing NumPy and SciPy.

## Repository guide

- [`src/`](../src/): current C++ solver, parser and CLI.
- [`benchmarks/`](../benchmarks/): reference checks and benchmark drivers; [`benchmarks/results/`](../benchmarks/results/) holds the completed CPU wave CSVs.
- [`results/`](../results/): retained C++ and historical ledgers. Historical files are not current-engine claims.
- [`gpu/`](../gpu/): CUDA kernels and approximate PDHG prototype.
- [`reproduction/`](../reproduction/): pinned C++ Netlib reproduction.
- [`examples/`](../examples/): synthetic refinery examples, without operational-data or savings validation.
- [`engine/`](../engine/), [`parsers/`](../parsers/): historical Python reference solver and MPS readers, separate from the current C++ engine; used by `examples/` and a few `benchmarks/` checks.
- Retired material (single-file fallback engine `cpp-engine/` with its milestone snapshots and logs, and the Python `milp/` and `qp/` prototypes) lives on the `archive/legacy-engines-2026-10-08` branch. Its results do not apply to `src/`. The MILP random-model generator moved to [`tools/milp_check.py`](../tools/milp_check.py) and its edge-case models to [`tests/milp_edge/`](../tests/milp_edge/).

## License and access

Proprietary. Only official competition judges and organizers may read this submission for evaluation. No copying, modification, redistribution, commercial use or AI/ML use; see [LICENSE](../LICENSE) for full terms. The repository remains private until its owner chooses otherwise.

## Synthetic refinery warm starts

[39 what-if LP rows](../results/refinery_warmstart_58f77af/warmstart_table.csv), using the evaluated source snapshot `58f77af`, record independent reference validation of both solves without and with a previous basis at 1e-6, with zero recorded wrong answers. Median warm/cold iteration ratio is 0.302, recomputed from the CSV. This is synthetic-data iteration evidence only: reuse of the previous basis is inferred from diagnostic output, and millisecond wall times do not support a stable speedup claim. [Protocol, summary and provenance](../results/refinery_warmstart_58f77af/).

## MIPLIB widening, October 3, 2026

[40-case CSV and provenance](../results/miplib_widening_20261003/): six objective
matches, seven infeasibility matches, 26 engine time limits and one no-JSON
anomaly, with zero rows labeled wrong under the protocol. Evaluated source snapshot `58f77af`;
requested soft caps were engine 120 seconds and HiGHS 300 seconds. The current
implementation times out on the 26 cases. HiGHS also hit its separate cap on
11 of those. The `neos-1425699` cause of the missing machine-readable result is unconfirmed. Do not add this
selection to older MIPLIB counts without checking overlap.


## Published production and logistics applications

Two additional textbook LP cases complement the existing refinery example:
[steel production planning and Dantzig transportation](../examples/literature_lp/).
The evaluated source snapshot `34ea8b4` returns optimal 192000 illustrative
dollars/week (3 iterations) and 153.675 thousand dollars (7 iterations).
Original-model row/bound violations and reference objective discrepancies are
zero in the recorded checks. Complete points, hashes, citations and numerical
limits are retained. These small LPs do not establish industrial savings,
MILP scheduling, power dispatch or complete application coverage.
