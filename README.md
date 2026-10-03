![TARAL-LP - A sovereign solver core for refinery planning.](docs/assets/TARAL-LP.png)

[Official website](https://taral-lp.vercel.app/) · [Canonical Kaggle notebook](https://www.kaggle.com/code/hexraei/notebookeb52168008)

[![C++17 core](docs/assets/badge-cpp.svg)](src/) [![Netlib measured 93/93](docs/assets/badge-netlib.svg)](results/netlib_postmerge_58f77af/) [![MIPLIB 40 measured](docs/assets/badge-miplib.svg)](results/miplib_widening_20261003/) [![QPLIB 19 measured](docs/assets/badge-qplib.svg)](results/qplib_widening_20261003/) [![Website](docs/assets/badge-website.svg)](https://taral-lp.vercel.app/) [![Kaggle notebook](docs/assets/badge-kaggle.svg)](https://www.kaggle.com/code/hexraei/notebookeb52168008) [![Proprietary license](docs/assets/badge-license.svg)](LICENSE)

TARAL-LP is a from-scratch C++17 solver core for sparse industrial optimization: CPU linear programming (LP), mixed-integer linear programming (MILP) and convex quadratic programming (QP), with an approximate CUDA primal-dual hybrid gradient (PDHG) method for large continuous LPs. The evaluated CPU engine uses only the C++ standard library; HiGHS, SciPy and NumPy belong to benchmark tools, not the optimization engine. Experimental, not a production replacement for established solvers.

## Requirement alignment

Checked means **measured in the linked scope**, not complete industrial coverage. Unchecked work is a prototype, planned or not started.

<details open>
<summary>17 requirements, each mapped to evidence or implementation</summary>

- [x] **R1 - LP, MILP and QP:** measured [Netlib](results/netlib_postmerge_58f77af/), [selected MIPLIB](benchmarks/results/miplib_ledger.csv) and [QP measurements](results/qplib_widening_20261003/); limits remain.
- [ ] **R2 - Future solver classes:** prototype [module boundaries](src/taral.hpp); mixed-integer quadratic, nonlinear and mixed-integer nonlinear programming not started.
- [x] **R3 - Continuous methods:** measured [simplex, dual and interior-point paths](results/adversarial_4517277/); historical source snapshot, not a universal pass.
- [ ] **R4 - Mixed-integer search:** prototype [branch and bound with propagation](src/milp.cpp); full cutting-plane, presolve and heuristic coverage planned.
- [ ] **R5 - Sparse algebra:** prototype [own sparse factorization](src/); [34-instance scaling profile](results/r5_factorization_20261003/r5_summary.md) measured, broader scalability work remains.
- [ ] **R6 - Multi-core solve:** planned; evaluated CPU engine is single-threaded.
- [x] **R7 - GPU method:** measured [approximate PDHG on T4](results/gpu_pdhg/t4_industrial_20261003/); no exact-solver speed advantage claimed.
- [x] **R8 - Numerical checks:** measured [Netlib checks against the original constraints](results/netlib_postmerge_58f77af/); broader reliable convergence remains unverified.
- [x] **R9 - Independent solver core:** [C++ dependency inspection, warning-free build and independent validation](docs/foundation.md); reference solvers are test-only, not proof of full correctness.
- [x] **R10 - Large models:** measured [GPU LPs up to 1,092,610 decision variables](results/gpu_pdhg/t4_industrial_20261003/ledger.csv), approximately; not consistent exact industrial-scale solves.
- [x] **R11 - Difficult numerical models:** measured [seeded adversarial LP/MILP](results/adversarial_4517277/tests2_base_full.csv); failures and unresolved cases retained.
- [x] **R12 - Basic interface:** [documented command-line and C++ interfaces](docs/interface.md), measured 31/31 command-line checks and one LP API example; stable versioning and broader integration coverage remain.
- [x] **R13 - Standard benchmarks:** measured [Netlib](results/netlib_postmerge_58f77af/) and [selected MIPLIB](benchmarks/results/miplib_ledger.csv); broader benchmark coverage planned.
- [x] **R14 - Established-solver comparison:** measured [HiGHS objective and feasibility checks](results/netlib_postmerge_58f77af/); no general speed claim.
- [x] **R15 - Challenging large models:** measured [historical large-model stress](results/adversarial_4517277/tests2_methods_dual_full.csv); weak-relaxation industrial coverage remains incomplete.
- [x] **R16 - Inspectable, extensible foundation:** [documented module layout, independent checks and measured API use](docs/foundation.md); proprietary evaluation access only, broader extension validation remains.
- [ ] **R17 - Representative applications:** measured [published production and transportation LPs](examples/literature_lp/) alongside [refinery examples](examples/); mixed-integer scheduling, dispatch and broader application evidence remain planned.

</details>

## What the evidence shows

- **CPU solution validation:** [93/93 Netlib passes](results/netlib_postmerge_58f77af/) using both primal and dual simplex; 92/93 meet the tighter check, with GREENBEA the only case outside the tighter tolerance.
- **Large approximate LPs:** [rail4284](results/gpu_pdhg/t4_industrial_20261003/ledger.csv), with 1,092,610 decision variables, reports an approximate answer (`near_optimal`) at 1e-4 in 143.4s on a shared T4. Not an exact optimum certificate.
- **Failures stay visible:** [historical per-case stress results](results/adversarial_4517277/) retain objective failures, cases with no verified answer and reference disagreements alongside passes.

## Benchmarks

| Measurement | Result | Evidence |
| --- | --- | --- |
| CPU Netlib, primal + dual | **93/93 each**, 92/93 at the tighter 1e-8 check, no wrong answers; evaluated source snapshot `58f77af` | [Primal](results/netlib_postmerge_58f77af/primal_ledger.csv), [dual](results/netlib_postmerge_58f77af/dual_ledger.csv) |
| T4 PDHG, rail2586 / rail4284 | Approximate answer (`near_optimal`) at stopping tolerance 1e-4: **110.4s / 143.4s**; both hit the 300s cap at 1e-6 | [Per-case results and test conditions](results/gpu_pdhg/t4_industrial_20261003/) |
| T4 PDHG, PDS-100 | Approximate answer (`near_optimal`) at stopping tolerance 1e-4 in **8.7s**, but objective error **7.72e-4** and violation of the original constraints **1.02** | [CSV](results/gpu_pdhg/t4_industrial_20261003/ledger.csv) |
| MIPLIB, 40 measured cases | **6 objective + 7 infeasibility matches**; 26 time limits, 1 run without a machine-readable result, zero rows labeled wrong | [CSV and limits](results/miplib_widening_20261003/) |
| QPLIB, 19 measured candidates | **5 objective matches + 1 locally verified better point; 6 feasible engine-reported optima** where HiGHS had no optimal reference; [investigated discrepancies/limits](docs/limitations.md#qplib-measurements) | [Per-case results and test conditions](results/qplib_widening_20261003/) |
| Historical adversarial LP, dual | **2903/2952 accepted**, including 29 accepted results with independently checked reference discrepancies; 15 objective failures, 34 cases with no verified answer; evaluated source snapshot `4517277` | [CSV](results/adversarial_4517277/tests2_methods_dual_full.csv) |

Protocol: Netlib engine time limit 60s/case, PILOT.WE and PILOT4 included with HiGHS references capped at 300s; default `reproduction/reproduce.sh` excludes `pilot.we`/`pilot4` because their reference answers are excluded in that protocol (`reference_excluded_60s`) (at most 91/93 counted passes).

The tighter CPU check requires row and bound violations and relative objective error at or below 1e-8. The ordinary pass rule is documented in [full measurements](docs/measured-results.md#verification-and-reproduction).

GPU results are approximate, double precision, one run on a shared T4. Solver tolerance is not a bound on objective error or violation of the original constraints; rail objective references are unavailable. [Full measurements and historical tables](docs/measured-results.md).

## Sparse algebra scaling profile

[Sparse LU measurements](results/r5_factorization_20261003/r5_factorization.csv) and [full summary](results/r5_factorization_20261003/r5_summary.md): dual simplex, 34 benchmark instances, 300 s cap, one run each, source commit 6215801. No singular or failed factorizations were measured. Fill ratio is about 1.0 on Kennington/OSA/PDS structure (maximum 1.04 on PDS), and refactorization time grows about linearly with basis size.

The large-model limit in this measurement is the per-iteration forward/backward solves: at 40,000 rows and above they take 21% to 74% of wall time. Wall time per iteration grows from 0.14 ms at 2,426 rows to 8 ms at 156,243 rows. The measured implementation does not yet exploit right-hand-side sparsity in those solves. This is a scaling profile, not a speed advantage.

21 instances solved to optimality and match HiGHS within 1e-6 relative; 13 reach the 300 s cap with no answer, while HiGHS finished all of them. PDS-30 refactorizes about three times as often as its neighbors; the cause was not isolated. These historical measurements do not close R5's broader scalability work.

## Reproduce in three commands

From the repository root, with Python 3, g++ and git installed:

```bash
python3 -m venv .venv && . .venv/bin/activate
python3 -m pip install -r reproduction/requirements.txt
bash reproduction/reproduce.sh
```

This builds the CPU engine, verifies the fixed benchmark files and runs the default 60s validation plus the additional TRUSS case, not the headline protocol that includes PILOT.WE and PILOT4. Counts can vary with host speed. [Reproduction details](reproduction/README_REPRO.md) · [Build and smoke checks](benchmarks/SMOKE.md) · [API](src/taral.hpp)

## Team

**Team IRIZ_**

<p>
<a href="https://github.com/Hexraei" title="Navin"><img src="docs/assets/team/Hexraei.png" width="56" height="56" alt="Navin" title="Navin" style="border-radius:50%;"></a>
<a href="https://github.com/muffedd" title="Sutharshan"><img src="docs/assets/team/muffedd.png" width="56" height="56" alt="Sutharshan" title="Sutharshan" style="border-radius:50%;"></a>
<a href="https://github.com/abinaya2006" title="Abinaya"><img src="docs/assets/team/abinaya2006.png" width="56" height="56" alt="Abinaya" title="Abinaya" style="border-radius:50%;"></a>
<a href="https://github.com/ksm-13" title="Manaswini"><img src="docs/assets/team/ksm-13.png" width="56" height="56" alt="Manaswini" title="Manaswini" style="border-radius:50%;"></a>
<a href="https://github.com/Professor-Mady" title="Madhesh"><img src="docs/assets/team/Professor-Mady.png" width="56" height="56" alt="Madhesh" title="Madhesh" style="border-radius:50%;"></a>
<a href="https://github.com/Saibts" title="Sailakshmi"><img src="docs/assets/team/Saibts.png" width="56" height="56" alt="Sailakshmi" title="Sailakshmi" style="border-radius:50%;"></a>
</p>

## Scope and access

[Limitations, unsupported models and numerical caveats](docs/limitations.md) · [Detailed results](docs/measured-results.md) · [License](LICENSE)

Proprietary submission. Access is restricted to official competition judges and organizers for evaluation. This README does not change repository visibility or grant reuse rights.
