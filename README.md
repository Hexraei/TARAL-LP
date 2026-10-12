![TARAL-LP - A sovereign solver core for refinery planning.](docs/assets/TARAL-LP.png)

[Official website](https://taral-lp.vercel.app/) · [Canonical Kaggle notebook](https://www.kaggle.com/code/hexraei/taral-lp-prototype-and-measured-benchmarks)

[![C++17 core](docs/assets/badge-cpp.svg)](src/) [![Netlib measured 93/93](docs/assets/badge-netlib.svg)](results/netlib_postmerge_58f77af/) [![MIPLIB 39/65 audited](docs/assets/badge-miplib.svg)](results/milp_paired_20261012/) [![QPLIB 19 measured](docs/assets/badge-qplib.svg)](results/qplib_widening_20261003/) [![Website](docs/assets/badge-website.svg)](https://taral-lp.vercel.app/) [![Kaggle notebook](docs/assets/badge-kaggle.svg)](https://www.kaggle.com/code/hexraei/taral-lp-prototype-and-measured-benchmarks) [![Proprietary license](docs/assets/badge-license.svg)](LICENSE)

TARAL-LP is a from-scratch C++17 solver core for sparse industrial optimization: CPU linear programming (LP), mixed-integer linear programming (MILP) and convex quadratic programming (QP), with an approximate CUDA primal-dual hybrid gradient (PDHG) method for large continuous LPs. The evaluated CPU engine uses only the C++ standard library; HiGHS, SciPy and NumPy belong to benchmark tools, not the optimization engine. Experimental, not a production replacement for established solvers.

## Measured results - October 12, 2026

**93/93 Netlib LP passes in each simplex mode. 39/65 MIPLIB 3 solver-optimal results in the latest independently audited experimental measurement.**

All 65 mixed-integer cases were measured in three paired repeats at 300 seconds per case. Adding objective-integrality bound strengthening raises the expected-value count from 38 to 39 in every repeat, with `vpm1` the sole new solve. All 366 returned points pass independent checks on the original models. [Per-case results and conditions](results/milp_paired_20261012/).

<sub>Numerical validation, not an exact optimality certificate; `vpm1` retains an inconclusive numerical cut check. At the 270-second margin the paired counts are 37/38; frozen-reference counts are 37/38. The 39/65 measurement is an experimental snapshot, not the unmodified main build.</sub>

## What sets TARAL-LP apart

- Own C++17 optimization engine: existing solvers are used for test references, not inside the solve path.
- Sparse primal and dual simplex, convex-QP interior-point method, and mixed-integer search, with original-model answer checks.
- Inspectable per-case measurements, separate unresolved outcomes, and an approximate GPU LP method with accuracy limits stated.

## Requirement alignment

Checked means measured or implemented in the named scope, not complete industrial coverage. Experimental measurements and the checked-in engine are distinguished below.

<details open>
<summary>17 requirements, each mapped to evidence or implementation</summary>

- [x] **LP, MILP and QP:** measured [Netlib](results/netlib_postmerge_58f77af/), [all 65 MIPLIB 3 cases](results/milp_paired_20261012/) and [QP measurements](results/qplib_widening_20261003/); limits remain.
- [x] **Future solver classes:** non-shipping [convex MIQP experiment](experiments/r2_miqp/) (diagonal convex Q, minimization); [original measured package](results/r2_miqp_20261008/SUMMARY.md): 1000/1000 seeded random and 32/32 structured cases agree with exhaustive-enumeration/SLSQP reference within 1e-6 relative (max 5.4e-07). Separate review fix prevents false optimal claims on unresolved subtrees; those paths were not covered by the original corpus. No presolve, warm starts or cuts; nonlinear/MINLP not started.
- [x] **Continuous methods:** measured [simplex, dual and interior-point paths](results/adversarial_4517277/); historical source snapshot, not a universal pass.
- [x] **Mixed-integer search:** [latest paired measurement](results/milp_paired_20261012/): 39/65 solver-optimal with original-model point checks and expected-value agreement in all three repeats; 38 finish within 270 seconds. Experimental cuts and objective-integrality strengthening are separate from the checked-in [branch-and-bound engine](src/milp.cpp). [Earlier 65-case baseline](results/r4_miplib3_20261008/) remains historical evidence.
- [x] **Sparse algebra:** prototype [own sparse factorization](src/); [34-instance scaling profile](results/r5_factorization_20261003/r5_summary.md) measured: fill ratio ~1.0 on network/PDS structure, mean refactorization about 1 microsecond per row, 21/34 solved optimal all matching HiGHS within 1e-6 relative, 13/34 timed out at the 300 s cap (HiGHS finished all of them); broader scalability work remains.
- [ ] **Multi-core solve:** planned; evaluated CPU engine is single-threaded.
- [x] **GPU method:** measured [approximate PDHG on T4](results/gpu_pdhg/t4_industrial_20261003/); no exact-solver speed advantage claimed.
- [x] **Numerical checks:** measured [Netlib checks against the original constraints](results/netlib_postmerge_58f77af/); broader reliable convergence remains unverified.
- [x] **Independent solver core:** [C++ dependency inspection, warning-free build and independent validation](docs/foundation.md); reference solvers are test-only, not proof of full correctness.
- [x] **Large models:** measured [GPU LPs up to 1,092,610 decision variables](results/gpu_pdhg/t4_industrial_20261003/ledger.csv), approximately; not consistent exact industrial-scale solves.
- [x] **Difficult numerical models:** measured [seeded adversarial LP/MILP](results/adversarial_4517277/tests2_base_full.csv); failures and unresolved cases retained.
- [x] **Basic interface:** [documented command-line and C++ interfaces](docs/interface.md), measured 31/31 command-line checks and one LP API example; stable versioning and broader integration coverage remain.
- [x] **Standard benchmarks:** measured [Netlib](results/netlib_postmerge_58f77af/) and [all 65 MIPLIB 3 cases](results/milp_paired_20261012/); broader benchmark coverage planned.
- [x] **Established-solver comparison:** measured [HiGHS objective and feasibility checks](results/netlib_postmerge_58f77af/); no general speed claim.
- [x] **Challenging large models:** measured [historical large-model stress](results/adversarial_4517277/tests2_methods_dual_full.csv); weak-relaxation industrial coverage remains incomplete.
- [x] **Inspectable, extensible foundation:** [documented module layout, independent checks and measured API use](docs/foundation.md); proprietary evaluation access only, broader extension validation remains.
- [ ] **Representative applications:** [published production and transportation LPs](examples/literature_lp/), refinery illustrations, and small [unit-commitment](examples/unit_commitment/) and [shift-scheduling](examples/shift_scheduling/) examples with solver-optimal status and reference objective agreement. These are illustrations, not industrial deployment evidence; the two MILP examples lack independently retained point checks.

</details>

## What the evidence shows

- **CPU solution validation:** [93/93 Netlib passes](results/netlib_postmerge_58f77af/) using both primal and dual simplex; 92/93 meet the tighter check, with GREENBEA the only case outside the tighter tolerance.
- **Large approximate LPs, not verified optimal solves:** [rail4284](results/gpu_pdhg/t4_industrial_20261003/ledger.csv), with 1,092,610 decision variables, stopped as `near_optimal` at tolerance 1e-4 in 143 s on a T4, with original-row violation **1.2e-3**. A later K=1 row-feasible run (violation under 1e-6) stopped as `near_optimal` in 2677 s on a T4. These figures come from the archived run-log transcript; the raw logs are no longer available. Objective about 1054.06 is PDHG's own value, not an independently verified optimum. No independent reference objective exists: HiGHS hit its 600 s cap, and the CPU simplex engine timed out at 300 s.
- **Failures stay visible:** [historical per-case stress results](results/adversarial_4517277/) retain objective failures, cases with no verified answer and reference disagreements alongside passes.

## Benchmarks

| Measurement | Result | Evidence |
| --- | --- | --- |
| CPU Netlib, primal + dual | **93/93 each**, 92/93 at the tighter 1e-8 check, no wrong answers; evaluated source snapshot `58f77af` | [Primal](results/netlib_postmerge_58f77af/primal_ledger.csv), [dual](results/netlib_postmerge_58f77af/dual_ledger.csv) |
| T4 PDHG, rail2586 / rail4284 | **Approximate, not verified optimal solves.** At stopping tolerance 1e-4: rail2586 110.4 s (original-row violation 9.7e-4); rail4284 143 s (original-row violation **1.2e-3**). Both original 1e-6 runs hit the 300 s cap. A later rail4284 K=1 row-feasible `near_optimal` run (violation under 1e-6) took 2677 s on a T4. Later-run figures are from the archived run-log transcript, not retained raw logs; no independent reference objective for rail4284. | [Per-case results and test conditions for the original runs](results/gpu_pdhg/t4_industrial_20261003/) |
| T4 PDHG, PDS-100 | Approximate answer (`near_optimal`) at stopping tolerance 1e-4 in **8.7s**, but objective error **7.72e-4** and violation of the original constraints **1.02** | [CSV](results/gpu_pdhg/t4_industrial_20261003/ledger.csv) |
| Latest MIPLIB 3, 65 cases x 3 paired repeats | **39/65 expected-value solves**, 38 within the 270s margin; all returned points independently checked. Experimental snapshot, not main. | [Ledger and numerical caveats](results/milp_paired_20261012/) |
| Historical MIPLIB selection, 40 measured cases | **6 objective + 7 infeasibility matches**; 26 time limits, 1 run without a machine-readable result, zero rows labeled wrong | [CSV and limits](results/miplib_widening_20261003/) |
| Historical QPLIB selection, 19 measured candidates | **5 objective matches + 1 locally verified better point; 6 feasible engine-reported optima** where HiGHS had no optimal reference; [investigated discrepancies/limits](docs/limitations.md#qplib-measurements) | [Per-case results and test conditions](results/qplib_widening_20261003/) |
| Historical adversarial LP, dual | **2903/2952 accepted**, including 29 accepted results with independently checked reference discrepancies; 15 objective failures, 34 cases with no verified answer; evaluated source snapshot `4517277` | [CSV](results/adversarial_4517277/tests2_methods_dual_full.csv) |

Protocol: Netlib engine time limit 60s/case, PILOT.WE and PILOT4 included with HiGHS references capped at 300s; default `reproduction/reproduce.sh` excludes `pilot.we`/`pilot4` because their reference answers are excluded in that protocol (`reference_excluded_60s`) (at most 91/93 counted passes).

The tighter CPU check requires row and bound violations and relative objective error at or below 1e-8. The ordinary pass rule is documented in [full measurements](docs/measured-results.md#verification-and-reproduction).

GPU results are approximate, double precision, one run on a shared T4. Solver tolerance is not a bound on objective error or violation of the original constraints; rail objective references are unavailable. [Full measurements and historical tables](docs/measured-results.md).

## Sparse algebra scaling profile

[Sparse LU measurements](results/r5_factorization_20261003/r5_factorization.csv) and [full summary](results/r5_factorization_20261003/r5_summary.md): dual simplex, 34 benchmark instances, 300 s cap, one run each, source commit 6215801. No singular or failed factorizations were measured. Fill ratio is about 1.0 on Kennington/OSA/PDS structure (maximum 1.04 on PDS), and refactorization time grows about linearly with basis size.

The large-model limit in this measurement is the per-iteration forward/backward solves: at 40,000 rows and above they take 21% to 74% of wall time. Wall time per iteration grows from 0.14 ms at 2,426 rows to 8 ms at 156,243 rows. The measured implementation does not yet exploit right-hand-side sparsity in those solves. This is a scaling profile, not a speed advantage.

21 instances solved to optimality and match HiGHS within 1e-6 relative; 13 reach the 300 s cap with no answer, while HiGHS finished all of them. PDS-30 refactorizes about three times as often as its neighbors; the cause was not isolated. These historical measurements do not establish broader industrial scalability.

## Reproduce the checked-in LP measurement in three commands

From the repository root, with Python 3, g++ and git installed:

```bash
python3 -m venv .venv && . .venv/bin/activate
python3 -m pip install -r reproduction/requirements.txt
bash reproduction/reproduce.sh
```

This does not reproduce the separate latest mixed-integer experimental snapshot. It builds the checked-in CPU engine, verifies the fixed benchmark files and runs the default 60s validation plus the additional TRUSS case, not the headline protocol that includes PILOT.WE and PILOT4. Counts can vary with host speed. [Reproduction details](reproduction/README_REPRO.md) · [Build and smoke checks](benchmarks/SMOKE.md) · [API](src/taral.hpp)

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

Proprietary submission. Evaluation rights are limited to official competition judges and organizers; the repository is publicly viewable. This README does not change repository visibility or grant reuse rights.
