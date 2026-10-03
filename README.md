![TARAL-LP - A sovereign solver core for refinery planning.](docs/assets/TARAL-LP.png)

[![C++17 core](docs/assets/badge-cpp.svg)](src/) [![Netlib 93/93](docs/assets/badge-netlib.svg)](results/netlib_postmerge_58f77af/) [![Proprietary license](docs/assets/badge-license.svg)](LICENSE)

TARAL-LP is a from-scratch C++17 solver core for sparse industrial optimization: CPU LP, MILP and convex QP, with an approximate CUDA PDHG path for large continuous LPs. The judged CPU engine uses only the C++ standard library; HiGHS, SciPy and NumPy belong to benchmark tools, not the solve path. Experimental, not a production replacement for established solvers.

## Requirement alignment

Checked means **measured in the linked scope**, not complete industrial coverage. Unchecked work is a prototype, planned or not started.

<details open>
<summary>17 requirements, each mapped to evidence or implementation</summary>

- [x] **R1 - LP, MILP and QP:** measured [Netlib](results/netlib_postmerge_58f77af/), [selected MIPLIB](benchmarks/results/miplib_ledger.csv) and [QP measurements](results/qplib_widening_20261003/); limits remain.
- [ ] **R2 - Future solver classes:** prototype [module boundaries](src/taral.hpp); MIQP, NLP and MINLP not started.
- [x] **R3 - Continuous methods:** measured [simplex, dual and interior-point paths](results/adversarial_4517277/); historical source snapshot, not a universal pass.
- [ ] **R4 - Mixed-integer search:** prototype [branch and bound with propagation](src/milp.cpp); full cutting-plane, presolve and heuristic coverage planned.
- [ ] **R5 - Sparse algebra:** prototype [own sparse factorization](src/); broader scalability validation planned.
- [ ] **R6 - Multi-core solve:** planned; headline CPU engine is single-threaded.
- [x] **R7 - GPU path:** measured [approximate PDHG on T4](results/gpu_pdhg/t4_industrial_20261003/); no exact-solver speed advantage claimed.
- [x] **R8 - Numerical checks:** measured [original-model Netlib gate](results/netlib_postmerge_58f77af/); broader reliable convergence remains unverified.
- [ ] **R9 - Independent solver core:** prototype [standard-library-only C++ engine](src/), [build definition](CMakeLists.txt); reference solvers are test-only.
- [x] **R10 - Large models:** measured [GPU LPs up to 1,092,610 columns](results/gpu_pdhg/t4_industrial_20261003/ledger.csv), approximately; not consistent exact industrial-scale solves.
- [x] **R11 - Difficult numerical models:** measured [seeded adversarial LP/MILP](results/adversarial_4517277/tests2_base_full.csv); failures and unresolved cases retained.
- [ ] **R12 - Basic interface:** prototype [CLI](src/main.cpp) and [C++ API](src/taral.hpp); no GUI required for use.
- [x] **R13 - Standard benchmarks:** measured [Netlib](results/netlib_postmerge_58f77af/) and [selected MIPLIB](benchmarks/results/miplib_ledger.csv); broader benchmark coverage planned.
- [x] **R14 - Established-solver comparison:** measured [HiGHS objective and feasibility checks](results/netlib_postmerge_58f77af/); no general speed claim.
- [x] **R15 - Challenging large models:** measured [historical large-model stress](results/adversarial_4517277/tests2_methods_dual_full.csv); weak-relaxation industrial coverage remains incomplete.
- [ ] **R16 - Inspectable, extensible foundation:** prototype [source](src/) and [independent checks](benchmarks/orig_check.py); see [access terms](LICENSE).
- [ ] **R17 - Representative applications:** prototype [refinery examples](examples/); [dataset scope](docs/limitations.md#qplib-scope) is documented, broader literature case studies planned.

</details>

## What the evidence shows

- **Exact CPU gate:** [93/93 Netlib passes](results/netlib_postmerge_58f77af/) in both modes; 92/93 meet the tighter check, with GREENBEA the only nonstrict case.
- **Large approximate LPs:** [rail4284](results/gpu_pdhg/t4_industrial_20261003/ledger.csv), with 1,092,610 columns, reports `near_optimal` at 1e-4 in 143.4s on a shared T4. Not an exact optimum certificate.
- **Failures stay visible:** [historical stress ledgers](results/adversarial_4517277/) retain objective failures, no-answer cases and reference disagreements alongside passes.

## Benchmarks

| Measurement | Result | Evidence |
| --- | --- | --- |
| CPU Netlib, primal + dual | **93/93 each**, 92/93 strict, no wrong answers; source `58f77af` | [Primal](results/netlib_postmerge_58f77af/primal_ledger.csv), [dual](results/netlib_postmerge_58f77af/dual_ledger.csv) |
| T4 PDHG, rail2586 / rail4284 | 1e-4 `near_optimal`: **110.4s / 143.4s**; both hit the 300s cap at 1e-6 | [CSV and protocol](results/gpu_pdhg/t4_industrial_20261003/) |
| T4 PDHG, PDS-100 | 1e-4 `near_optimal` in **8.7s**, but objective error **7.72e-4** and original-row residual **1.02** | [CSV](results/gpu_pdhg/t4_industrial_20261003/ledger.csv) |
| QPLIB, 19 measured candidates | **5 objective matches + 1 locally verified better point; 6 feasible engine-reported optima** where HiGHS had no optimal reference; [investigated discrepancies/limits](docs/limitations.md#qplib-measurements) | [CSV and protocol](results/qplib_widening_20261003/) |
| Historical adversarial LP, dual | **2903/2952 accepted**, including 29 reference-corrected passes; 15 objective failures, 34 no-answer cases; source `4517277` | [CSV](results/adversarial_4517277/tests2_methods_dual_full.csv) |

Protocol: Netlib engine cap 60s/case, pilots counted with HiGHS references capped at 300s; default `reproduction/reproduce.sh` excludes `pilot.we`/`pilot4` as `reference_excluded_60s` (ceiling 91/93).

GPU results are approximate, fp64, one run on a shared T4. Solver tolerance is not a bound on objective error or original-row residuals; rail objective references are unavailable. [Full measurements and historical tables](docs/measured-results.md).

## Reproduce in three commands

From the repository root, with Python 3, g++ and git installed:

```bash
python3 -m venv .venv && . .venv/bin/activate
python3 -m pip install -r reproduction/requirements.txt
bash reproduction/reproduce.sh
```

This builds the CPU engine, verifies the pinned corpus and runs the default 60s gate plus TRUSS, not the pilots-counted headline protocol. Counts can vary with host speed. [Reproduction details](reproduction/README_REPRO.md) · [Build and smoke checks](benchmarks/SMOKE.md) · [API](src/taral.hpp)

## Team

**Team IRIZ_**

| Member | GitHub |
| --- | --- |
| Navin | [Hexraei](https://github.com/Hexraei) |
| Sutharshan | [muffedd](https://github.com/muffedd) |
| Sailakshmi | [Saibts](https://github.com/Saibts) |
| Madhesh | [Professor-Mady](https://github.com/Professor-Mady) |
| Manaswini | [ksm-13](https://github.com/ksm-13) |
| Abinaya | [abinaya2006](https://github.com/abinaya2006) |

## Scope and access

[Limitations, unsupported models and numerical caveats](docs/limitations.md) · [Detailed results](docs/measured-results.md) · [License](LICENSE)

Proprietary submission. Access is restricted to official competition judges and organizers for evaluation. This README does not change repository visibility or grant reuse rights.
