# TARAL-LP

C++17 linear-programming solver with an MPS parser, sparse LU and bounded-variable revised simplex. The engine in `src/` uses the C++ standard library. HiGHS supplies benchmark reference results.

## Current status - October 2, 2026

Experimental CPU solver. Not ready for production; Kaggle verification is pending. Local `src/` measurements use engine commit `a9e8218` (primal ledger) and `0bad060` (dual ledger; LP verdicts of the default method are unchanged since) and 93 Netlib cases from the [pinned corpus](https://github.com/ozy4dm/lp-data-netlib/tree/56257eea85b433ce6aa67d26156b36385318fd6f/mps_files).

| Protocol | Result | Verification |
| --- | --- | --- |
| 60 seconds per case | 90/93 passes, 0 wrong answers | Local; no lost October 1 Python baseline passes |
| Stricter numerical check | 89/90 passing cases | Local; row/bound violation and relative objective error <= 1e-8 |
| 60 seconds, `--method dual` (separate line) | 91/93 passes, 0 wrong answers, STRICT 90/91 | Local; engine commit `0bad060`; DFL001 passes (about 14 s here); no case lost against the primal |
| 300 seconds per case | 92/93 passes, 0 wrong answers | Local; separate extended-cap protocol |
| MILP, 26 small MIPLIB instances | 19/26 solved at 300 s (14/26 at 30 s); 7 stop at the time limit, none wrong | 300 s run on a cloud container, 30 s run local; `results/milp_miplib/` |
| Convex QP (interior point), 99 Maros-Meszaros instances | 71 pass, 28 fail under the stated rule: 5 objective mismatches (checked: 4 are the HiGHS default reference being the less accurate side or misreading the file, the fifth is uncertified; see `results/qp_maros_meszaros/NOTES.md`), 18 without a HiGHS reference, 4 numerical failures, 1 nonconvex | Cloud container, 60 s; ledger rows unchanged |

At 60 seconds, DFL001 times out. PILOT.WE and PILOT4 lack usable HiGHS reference results in that protocol and remain excluded from the pass count; its ceiling is 91/93. They pass under the separate 300-second protocol; DFL001 still does not finish. TRUSS passes as an extra case outside the 93-case denominator.

A local check of 3,600 small LPs with known optimal/infeasible/unbounded status found 0 false verdicts.

Scope:
- LPs on CPU: primal simplex (default) and an opt-in dual simplex (`--method dual`). Branch-and-bound MILP, an interior-point QP path and `OBJSENSE` are implemented; see `results/milp_miplib`, `results/ipm_netlib`, `results/qp_maros_meszaros`.
- Dual line ledger: `results/cpp_0bad060_dual_60s`. It is reported beside, never merged with, the primal headline.
- Industrial-scale reliability remains untested.
- No measured speed advantage over HiGHS.
- GPU PDHG is a near-optimal prototype with its own ledgers in `results/gpu_pdhg`; Mittelmann and QPLIB are outside the tested coverage.

## Build and run

Build with a C++17 compiler. The engine has no external solver or linear-algebra dependency.

```bash
g++ -O3 -march=native -std=c++17 -o taral src/*.cpp
./taral path/to/model.mps --time-limit 60 --sol out.sol --json out.json
```

`out.json` reports status, objective including the objective-row RHS constant, iterations, wall time and a message. `out.sol` lists original variable values when status is optimal. Possible statuses are `optimal`, `infeasible`, `unbounded`, `time_limit`, `iteration_limit`, `numerical_failure` and `parse_error`.

## How results are checked

Fixed pass rule:
- Both TARAL-LP and original-MPS HiGHS return optimal status.
- Objective error <= `max(1e-6, 1e-7 * abs(HiGHS objective))`.
- Independently recomputed original objective agrees.
- Relative row violation <= `1e-6`, divided by `1 + abs(row RHS)`.
- Bound violation <= `1e-6`.

The stricter check is reported separately; the fixed rule stays unchanged.

The one 60-second pass outside the stricter line is GREENBEA: row `BRG...U3` has zero RHS and about 95 terms reaching 2.3e8 that cancel; measured violation is 1.46e-8. GREENBEA passes the fixed rule but fails the stricter check.

Local machine: Intel Core 7 240H, 10 cores/16 threads, 15 GiB RAM, Ubuntu 24.04, Linux 6.17, g++ 13.3 with `-O3 -march=native -std=c++17`. Cases ran one at a time on a single-threaded engine. Host speed affects near-cap outcomes. Use separate timing records for local and Kaggle runs. Kaggle verification of this `src/` engine is pending.

## Repository guide

| Path | Purpose |
| --- | --- |
| `src/` | Headline C++ engine: MPS parser, sparse LU, primal and dual simplex, branch-and-bound MILP, interior-point QP and CLI |
| `cpp-engine/` | Separate bundled C++ reference engine, diagnostic tools, fixtures and internal report; its results apply to this engine alone |
| `reproduction/` | One-command reproduction of the `src/` headline engine (build, pinned corpus, 93-case gate); see its [instructions](reproduction/README_REPRO.md) |
| `gpu/` | CUDA sparse matrix-vector kernel benchmark and a PDHG LP solver prototype (near-optimal, reported on separate ledgers in `results/gpu_pdhg/`); the exact LP headline is CPU-only |
| `benchmarks/` | Current benchmark and checking tooling for the C++ engine (`netlib_gate.py`, independent checker, MILP/QP/MPS-semantics suites; HiGHS is used only here) |
| `engine/`, `parsers/` | Historical Python/SciPy prototype; the current C++ engine is in `src/` |
| `milp/`, `qp/` | Historical small synthetic Python prototypes; the C++ MILP and QP results are the ledgers above |
| `results/` | Older partial Python ledgers; current C++ results are summarized above |
| `examples/` | Synthetic refinery LP and retained hypothetical objective; operational data and measured savings are unavailable |

Reproduction builds `src/` and runs `benchmarks/netlib_gate.py`; it requires a clone of this whole repository and fetches the pinned corpus. Keep `cpp-engine/` and historical results separate from the current `src/` headline.

## Historical measurements

<details>
<summary>Python prototype benchmarks and limitations (September 30 - October 1, 2026)</summary>

## Historical baseline - Python/SciPy prototype

The earlier Python prototype uses SciPy LU. These records document that prototype, outside the compliant C++ engine path. Measurements ran on Kaggle unless marked local.

- Ratios above 1 mean TARAL-LP took longer than HiGHS.
- The medians use different passing pairs. A controlled speed comparison requires identical pairs.
- The historical headline is the October 1 Python run.
- Four executed Kaggle runs of the same October 1 Python engine code gave 73, 73, 73 and 74 passes out of 93.
- The 73 October 1 Python passes held in every run, with no losses against the September 30 baseline.
- Host speed can change outcomes near the 60-second cap, including 80BAU3B.
- Median TARAL-LP / HiGHS solver-wall ratios ranged from about 4.67x to 5.93x, varying by up to about 25% between runs.
- Headline selection: the October 1 Python result, 73/93 and 5.828x. The best of the four runs is excluded from that headline.
- Both solvers had a 60-second parse + solve budget on Kaggle CPU.
- The Kaggle comparison uses the same parser for both solvers.
- Local original-MPS diagnostics checked all 73 passing cases of the same October 1 Python engine. This gate was limited to local passers; Kaggle and the 20 non-passing cases remained unchecked by it.
- LP solving remained CPU-only.

### Local benchmarks - patched Python engine (September 30, 2026, historical)

The patched Python engine kept 42/42 selected baseline cases passing and recovered 10 of the 35 cases that timed out in the earlier 57-pass Kaggle run. Scope: selected local checks only. The ten recoveries have no Kaggle confirmation; the full 93-case corpus was not checked locally.

Two changes were tested:

1. Faster candidate checks. Each simplex iteration builds a set of the variables already in the basis. Checking whether a candidate is in that set replaces repeated searches through a list. Candidate order, the reduced-cost threshold and tie-breaking stay unchanged.
2. Build the model matrix once. The parser lays out original rows, range rows and finite-bound rows before allocating and filling the final dense matrix. This removes repeated dictionary scans and whole-matrix copies. Row order, bound substitution and the objective offset stay unchanged. The parser output and solver are still dense; this does not remove large-model memory pressure.

#### Local test conditions

Intel Xeon @ 2.60 GHz, x86_64; two logical CPUs exposed by the OS (not a verified count of physical host cores); about 1.94 GiB RAM, no swap. Ubuntu 22.04.5 LTS, Linux 6.1.158+, Python 3.10.12, NumPy 2.2.6 and SciPy 1.15.3. Workers used `OPENBLAS_NUM_THREADS=1` and `OMP_NUM_THREADS=1`.

The recovery checks used up to six concurrent workers. Timing scope: concurrent-load wall time, including parsing and solving. Isolated CPU performance was not measured. The local per-case cap was 65 seconds; Kaggle's benchmark cap was 60 seconds. In particular, 25FV47 passed locally in 60.48 seconds, already outside the Kaggle budget.

- The 42 baseline regression cases passed against the saved HiGHS objective references.
- The ten recovered cases also matched their saved references within the test tolerance.
- The local ledger contains 77 checks: 42 baseline cases plus all 35 former timeouts.
- Untested by this ledger: the other 15 passes from the earlier 57-pass run and DEGEN2. Report the 52 local passes within their selected set; 52/93 would overstate coverage.

#### Ten former timeouts that passed locally

| Case | Local parse + solve (seconds) | Primal residual |
| --- | ---: | ---: |
| 25FV47 | 60.48 | 4.55e-13 |
| CZPROB | 11.77 | 9.09e-13 |
| FIT1P | 27.65 | 2.84e-14 |
| GANGES | 12.32 | 1.62e-11 |
| SCTAP2 | 35.13 | 2.33e-14 |
| SCTAP3 | 46.86 | 3.55e-15 |
| SHIP08L | 10.69 | 1.34e-14 |
| SHIP12L | 10.84 | 4.6e-14 |
| SIERRA | 55.93 | 4.55e-13 |
| STOCFOR2 | 39.92 | 4.55e-13 |

Residual: primal-constraint error under this test's calculation. Full optimality certification requires further checks.

#### What did not pass locally

Six solver failures:

- WOOD1P: iteration limit.
- MODSZK1: lost primal feasibility.
- PEROLD: iteration limit.
- MAROS: lost primal feasibility.
- PILOT4: lost primal feasibility.
- PILOTNOV: lost primal feasibility.

Fifteen 65-second caps: NESM, MAROS-R7, CYCLE, WOODW, PILOT.JA, D6CUBE, DEGEN3, BNL1, PILOT87, BNL2, PILOT.WE, PILOT, D2Q06C, GREENBEA, GREENBEB.

Four crashed/unclassified jobs: FIT2P, 80BAU3B, DFL001, FIT2D.

Crash causes remain unknown; both numerical-failure and recovery classifications remain unverified. DEGEN2 remains unresolved: the earlier 57-pass engine lost primal feasibility. The two speed fixes leave that issue unverified.

#### Why local times are not Kaggle times

A separate sequential, single-BLAS-thread check used the unchanged earlier 57-pass engine and the same pinned MPS files for 16 passing cases. Kaggle/local solver-wall time had a 1.974x median: Kaggle took about twice as long in this sample. The 10th-90th percentile range was 1.763-2.507x; the full range was 1.711-11.062x, with tiny AFIRO the overhead-heavy outlier.

- Environment comparison: hardware, runtime and load were not isolated.
- Runtime/library versions, startup and load can contribute.
- The timing comparison uses the unchanged earlier engine on both hosts.
- Hypothetical scaling of the concurrent-load recovery times, still awaiting an executed test: CZPROB, GANGES, SHIP08L and SHIP12L have room under 60 seconds; FIT1P is borderline; SCTAP2, STOCFOR2, SCTAP3, SIERRA and 25FV47 are at risk.
- Scaling parse time with a solver-time factor adds uncertainty.
- Use the saved September 30 run for the historical count. The scaling scenario carries forecast uncertainty.

### Executed Kaggle CPU benchmarks (Python prototype, historical)

Saved CPU runs; the notebook is private. HiGHS ran separately as the reference solver. Ratios above 1 mean TARAL-LP took longer.

#### Selected 42-case baseline

- 42/42 selected Netlib cases passed. Coverage: selected cases only.
- Median TARAL-LP/HiGHS wall-time ratio: 90.6x on the selected passing pairs.
- Maximum reported primal residual: 4.66e-10 across those passes.
- 30/30 small synthetic MILP tests and 30/30 small synthetic QP tests executed and passed. The MILP prototype is limited branch-and-bound; the QP prototype handles positive-definite convex box-bounded problems. Coverage excludes MIPLIB and general QP benchmarks.
- KKT optimality conditions were checked locally for 12/42 selected LPs. Certification of all 42 cases and Kaggle certificate execution remain pending.

#### Full-corpus run - October 1, 2026

- Attempted: 93 Netlib cases, pinned corpus commit `56257eea85b433ce6aa67d26156b36385318fd6f`.
- Result: 73 of 93 attempted cases passed, 13 timeouts, 5 `solver_failed`, 2 `highs_failed`.
- Median TARAL-LP/HiGHS solver-wall ratio: 5.828x over 73 passing pairs only. TARAL-LP is slower by this median factor. Earlier medians used different pair sets; a like-for-like speed comparison requires identical pairs.
- Maximum passing primal residual: 9.31e-10 (LOTFI).
- Saved-run ID: `354282101`. Figures in this section are from the October 1 Python result; the four-run host-speed variation check is described above. Both solvers had 60 seconds per case for parse + solve; paired ratios use solver-wall fields.
- Kaggle verification covers shared-parser formulation/objective match. Original-MPS diagnostics have a separate local gate, described below.

Named October 1 timeouts: 80BAU3B, D2Q06C, D6CUBE, DEGEN3, DFL001, FIT2D, FIT2P, GREENBEA, GREENBEB, PEROLD, PILOT, PILOT87, WOODW.

Named October 1 `solver_failed` cases: CYCLE, DEGEN2, MAROS, PILOT.JA, PILOTNOV.

- Named October 1 `highs_failed` cases: PILOT.WE and PILOT4.
- PEROLD is classified as timeout, but its HiGHS attempt also failed.
- Pass count contribution from these cases: zero.
- Three cases lack a usable HiGHS optimum in this run, leaving 90 of 93 comparable with HiGHS until reference/model issues are resolved.
- Headline denominator: 93 attempted. The 90 comparable cases and 73 selected successes are subsets.

#### Separate local original-MPS diagnostic - October 1 passers

- All 73 of 73 October 1 Python passing cases were independently checked locally using the same October 1 Python engine code. `benchmarks/orig_check.py` uses its own MPS parser, maps the returned solution back to original variables, checks row activities and bounds, and compares the recomputed original objective with HiGHS reading the original MPS file.
- Coverage: local only, limited to those passers.
- Kaggle and the 20 non-passing cases remain outside this diagnostic gate.

- The per-case pass rule requires HiGHS optimal status, objective error no greater than `max(1e-6, 1e-7 * abs(HiGHS objective))`, relative original-row violation no greater than `1e-6`, and absolute bound violation no greater than `1e-6`.
- Relative row violation is normalized by `1 + abs(original row RHS)`.
- The checker applies a 60-second solver limit per case.
- Observed maxima across the 73 cases: 1.05e-9 relative original-row violation and 9.3e-10 absolute bound violation.

- E226 objective convention: the MPS objective-row RHS entry is -7.113, which contributes an additive objective offset of +7.113.
- The solver-reported objective omits this offset: approximately -18.7519 versus the checker's and original-MPS HiGHS objective -11.6389.
- The corpus/published-reference convention excludes that offset, so E226 remains a pass under that convention; its solution vector is optimal and feasible in the independent check.
- The Python solver's reported objective therefore omits an original-MPS constant in this case.

To run the checker from the repository root with an existing MPS corpus and `highspy` installed:

```bash
python benchmarks/orig_check.py --engine engine/revised_simplex.py --corpus /path/to/mps_files --out /tmp/orig_check_out.json afiro e226
```

#### What changed from the September 30 baseline

- All 68 September 30 passing cases remained passing.
- The five gains are BNL1, BNL2, MAROS-R7, MODSZK1 and WOOD1P.
- Relative to the earlier 72/93 result, MAROS-R7 is newly passing.
- The 42/42 regression, 30/30 synthetic MILP and 30/30 synthetic QP checks remained green, with QP maximum objective error 3.197e-13.
- Separate selected gates are EXTRA_NETLIB 7/7 and CPU_LP 6/7; CPU_LP is not an all-pass result.
- The 20 named non-passing full-corpus outcomes remain.

- Earlier 72-pass run: 72/93 passing, 14 timeout, 5 solver_failed, 2 highs_failed, median 8.4156x over 72 passing pairs, maximum passing residual 8.15e-10.
- September 30 baseline: 68/93 passing, 16 timeout, 7 solver_failed, 2 highs_failed, median 12.478x over 68 pairs, maximum passing residual 9.31e-10.
- Earlier 57-pass run: 57/93 passing, 35 timeout, 1 solver_failed, median 29.98x over 57 pairs, maximum passing residual 9.31e-10.
- The passing sets differ, preventing a controlled speedup comparison.

- The full-corpus input source is [the pinned Netlib MPS collection](https://github.com/ozy4dm/lp-data-netlib/tree/56257eea85b433ce6aa67d26156b36385318fd6f/mps_files).
- The new local and Kaggle ledgers are separate from the older checked-in partial CSVs below.
- The checked-in engine is synced to the executed October 1 Python solve path, and the parser retains the September 30 bulk-construction path.
- Older checked-in partial ledgers are listed below; the new full-run ledger is separate.


## Older checked-in measurements (historical Python prototype)

- The consolidated partial sweep ledger `results/consolidated_partial_sweep.csv` records 39 distinct passing LP instances and five failed rows.
- A separate verified SCSD8 addendum (`results/scsd8_verified_addendum.jsonl`) raises that historical distinct verified count to 40.
- Scope: older retained measurements, separate from the local and Kaggle totals above.
- SHIP04S hit a run limit and remains excluded from verified passes.
- The retained results match HiGHS objectives and published Netlib reference objectives, with reported maximum primal residual 3.79e-10 across those passes.
- The ledger includes five failed rows; the full attempt inventory is incomplete.
- An independent rerun from the checked-in files remains pending.
- Reproduction of these older results requires missing handoff artifacts: original MPS inputs, reference source links, environment capture and a complete run script.
- End-to-end reproduction of this older sweep requires those artifacts.

FIT1D is a separate result file with objective -9146.378092420928 versus HiGHS -9146.378092420926 and published reference -9146.3780924; it took 42.45 CPU seconds versus HiGHS 0.0223 seconds.

Unverified attempts in that older sweep included SCORPION (singular basis), FINNIS and SCRS8 (parsed HiGHS result mismatches published references), SCAGR25 and MODSZK1 (lost primal feasibility), and reported outside these attached ledgers: DEGEN2, TUFF, SCTAP2, WOOD1P. Use the newer ledgers for later outcomes.

- The separate 14-case same-machine CSV applies only to the earlier pre-LU engine. Current LU performance needs its own measurements.
- It reports HiGHS faster on all 14 tested cases; the median ratio of our CPU time to HiGHS CPU time across the 14 rows is about 52.7x.
- This prototype was slower than HiGHS on the tested set.
- Timing scope: process CPU measurements on the supplied CSV cases. Performance elsewhere remains untested.
- HiGHS ran separately as the reference comparator.

- Historical LP solves ran on CPU; QP support is limited to the box-bounded prototype below.
- A separate small, pure-integer/binary MILP branch-and-bound prototype is included in `milp/`; its three fixed synthetic cases match HiGHS.
- A seeded property test adds 30/30 randomized small synthetic cases matching brute force and SciPy/HiGHS.
- Large-instance MILP reliability remains untested; small-MIPLIB results are in `results/milp_miplib/`.
- The separate QP prototype covers only positive-definite convex box-bounded problems; its three fixed synthetic cases and 30/30 seeded randomized SPD box-QP property tests match SciPy objectives.
- General QP and public benchmark coverage remain untested.
- Any future GPU matrix-operation experiments remain separate from LP solves.
- The CPU implementation uses sparse basis LU with product-form updates, plus dense fallback.
- Large refinery models remain out of present verified scope.
- Refinery examples are synthetic; the repository contains no refinery field data.



</details>

## Illustrative example

`examples/illustrative_refinery.py` is a synthetic six-variable LP. Its retained result matches HiGHS at about $1.19 million/day in assumed units. That hypothetical objective has no operational-data, measured-margin or savings validation. It uses the historical Python engine. Run from the repository root after installing NumPy and SciPy:

```bash
python -m examples.illustrative_refinery
```

## License and access

Proprietary. Only official competition judges and organizers may read this submission for evaluation. No copying, modification, redistribution, commercial use or AI/ML use; see [LICENSE](LICENSE) for full terms. The repository remains private until its owner chooses otherwise.
