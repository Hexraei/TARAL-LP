# TARAL-LP

TARAL-LP is proprietary. Only official competition judges and organizers may read it, solely to evaluate this submission.
No copying, changes, redistribution, commercial use, or use as AI/ML training data or input. No automated scraping or analysis except non-AI evaluation by competition organizers.
See LICENSE for the full terms.

An experimental CPU linear-programming prototype for refinery planning. It uses a two-phase revised simplex solver and reads MPS model files. It is not a production refinery planner or a replacement for HiGHS, CPLEX, or Xpress.

## Improved since Sep 30

Latest executed CPU result: **Kaggle v23, October 1, 2026** (saved-run ID `354282101`). This block is replaced when verified results change, not extended with daily logs.

| Measurement | Sep 30 baseline (v19) | Latest (v23) |
| --- | ---: | ---: |
| Passing cases / attempted | 68 / 93 | 73 / 93 |
| Timeouts | 16 | 13 |
| Solver-failed classifications | 7 | 5 |
| HiGHS-failed classifications | 2 | 2 |
| Maximum passing primal residual | 9.31e-10 | 9.31e-10 |
| Median TARAL-LP / HiGHS solver-wall ratio | 12.478x over 68 pairs | 5.828x over 73 pairs |

**Five added passes, no lost v19 baseline passes:** BNL1, BNL2 and MAROS-R7 changed from timeout to pass; MODSZK1 and WOOD1P changed from solver failure to pass. MAROS-R7 is the added pass since v21's 72/93. The selected regression remains **42/42**; executed small synthetic checks remain **30/30 MILP and 30/30 QP** (QP maximum objective error **3.197e-13**). These are limited prototypes, not general MILP/QP coverage. Other executed selected gates: **EXTRA_NETLIB 7/7**, **CPU_LP 6/7**. The latter is not all-green; these separate gates do not change the 93-case denominator.

Ratios above 1 mean TARAL-LP took longer than HiGHS. The paired sets changed, so the medians are **not a like-for-like speed comparison**. This is one executed run, not a repeated timing study. Both solvers had a 60-second parse + solve budget on Kaggle CPU. The Kaggle comparison uses the same parser for both solvers. A separate local original-MPS diagnostic run verified the same v23 engine on all 73 passing cases; it is not a Kaggle-executed or full-corpus verification gate and does not cover the 20 non-passing cases. No GPU LP solve is claimed.

## Local benchmarks - patched CPU engine (September 30, 2026)

A linear-programming solver finds the best value of an objective while respecting linear constraints. TARAL-LP is a prototype of that solver. The newest local checks kept **42/42 selected baseline cases passing** and recovered **10 of the 35 cases that timed out in Kaggle v18**. These are measured local results, not ten new Kaggle passes and not a full 93-case local result.

Two changes were tested:

1. **Faster candidate checks.** Each simplex iteration builds a set of the variables already in the basis. Checking whether a candidate is in that set replaces repeated searches through a list. Candidate order, the reduced-cost threshold and tie-breaking stay unchanged.
2. **Build the model matrix once.** The parser lays out original rows, range rows and finite-bound rows before allocating and filling the final dense matrix. This removes repeated dictionary scans and whole-matrix copies. Row order, bound substitution and the objective offset stay unchanged. The parser output and solver are still dense; this does not remove large-model memory pressure.

### Local test conditions

Intel Xeon @ 2.60 GHz, x86_64; two logical CPUs exposed by the OS (not a verified count of physical host cores); about 1.94 GiB RAM, no swap. Ubuntu 22.04.5 LTS, Linux 6.1.158+, Python 3.10.12, NumPy 2.2.6 and SciPy 1.15.3. Workers used `OPENBLAS_NUM_THREADS=1` and `OMP_NUM_THREADS=1`.

The recovery checks used up to six concurrent workers. The times below are **wall times under concurrent load, including parsing and solving**, not isolated CPU benchmarks. The local per-case cap was 65 seconds; Kaggle's benchmark cap was 60 seconds. In particular, 25FV47 passed locally in 60.48 seconds, already outside the Kaggle budget.

The 42 baseline regression cases passed against the saved HiGHS objective references. The ten recovered cases also matched their saved references within the test tolerance. The local ledger contains 77 checks: 42 baseline cases plus all 35 former timeouts. It does not retest the other 15 v18 passes or DEGEN2, so **52 local passes must not be reported as 52/93**.

### Ten former timeouts that passed locally

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

The residual measures how far the returned solution misses the primal constraints under this test's calculation. It is not a full optimality certificate.

### What did not pass locally

**Six solver failures:**

- WOOD1P: iteration limit.
- MODSZK1: lost primal feasibility.
- PEROLD: iteration limit.
- MAROS: lost primal feasibility.
- PILOT4: lost primal feasibility.
- PILOTNOV: lost primal feasibility.

**Fifteen 65-second caps:** NESM, MAROS-R7, CYCLE, WOODW, PILOT.JA, D6CUBE, DEGEN3, BNL1, PILOT87, BNL2, PILOT.WE, PILOT, D2Q06C, GREENBEA, GREENBEB.

**Four crashed/unclassified jobs:** FIT2P, 80BAU3B, DFL001, FIT2D.

The crashed jobs are not claimed as numerical failures or recoveries: their cause was not established by this ledger. **DEGEN2 remains a known issue:** v18 returned a primal-feasibility loss, and these two speed fixes do not establish a fix for it.

### Why local times are not Kaggle times

A separate sequential, single-BLAS-thread check used the unchanged v18 engine and the same pinned MPS files for 16 passing cases. Kaggle/local solver-wall time had a **1.974x median**: Kaggle took about twice as long in this sample. The 10th-90th percentile range was 1.763-2.507x; the full range was 1.711-11.062x, with tiny AFIRO the overhead-heavy outlier.

This is an environment comparison, not a hardware-only speed test. Runtime/library versions, startup and load can contribute. It does not mix the patched local engine with old Kaggle timings. Applying that factor to the concurrent-load recovery times is only a rough scenario, not an executed result: CZPROB, GANGES, SHIP08L and SHIP12L have room under 60 seconds; FIT1P is borderline; SCTAP2, STOCFOR2, SCTAP3, SIERRA and 25FV47 are at risk. Scaling parse time with a solver-time factor adds uncertainty. The saved v19 run settled that historical count; the scenario was not a guaranteed forecast.

## Executed Kaggle CPU benchmarks

These results come from saved, executed CPU runs. The notebook is private, so no notebook URL is published here. HiGHS is a separate reference solver, not part of the TARAL-LP solve path. A ratio above 1 means TARAL-LP took longer than HiGHS.

### v16 - selected 42-case baseline

- **42/42 selected Netlib cases passed.** This is a selected set, not the full feasible corpus.
- Median TARAL-LP/HiGHS wall-time ratio: **90.6x** on the selected passing pairs.
- Maximum reported primal residual: **4.66e-10** across those passes.
- **30/30 small synthetic MILP tests and 30/30 small synthetic QP tests executed and passed.** The MILP prototype is limited branch-and-bound; the QP prototype handles positive-definite convex box-bounded problems. These are not MIPLIB or general QP coverage.
- **12/42 local KKT checks** cover a subset of the selected LP set. KKT checks test optimality conditions; this local subset is not a 42-case certificate claim or an executed Kaggle certificate result.

### Latest completed full-corpus run - v23

- Attempted: **93 Netlib cases**, pinned corpus commit `56257eea85b433ce6aa67d26156b36385318fd6f`.
- Result: **73 of 93 attempted cases passed**, **13 timeouts**, **5 `solver_failed`**, **2 `highs_failed`**.
- Median TARAL-LP/HiGHS solver-wall ratio: **5.828x over 73 passing pairs only**. TARAL-LP is slower by this median factor. No like-for-like speed change is claimed against earlier medians with different pair sets.
- Maximum passing primal residual: **9.31e-10 (LOTFI)**.
- Saved-run ID: `354282101`. One executed run, no repeats yet. Both solvers had **60 seconds per case for parse + solve**; paired ratios use solver-wall fields.
- The Kaggle run checks shared-parser formulation/objective match. A separate local original-MPS diagnostic run is described below; it does not widen the Kaggle gate.

**Named v23 timeouts:** 80BAU3B, D2Q06C, D6CUBE, DEGEN3, DFL001, FIT2D, FIT2P, GREENBEA, GREENBEB, PEROLD, PILOT, PILOT87, WOODW.

**Named v23 `solver_failed` cases:** CYCLE, DEGEN2, MAROS, PILOT.JA, PILOTNOV.

**Named v23 `highs_failed` cases:** PILOT.WE and PILOT4. PEROLD is classified as timeout, but its HiGHS attempt also failed. None of these counts as a pass. Three cases lack a usable HiGHS optimum in this run, leaving **90 of 93 comparable with HiGHS** until reference/model issues are resolved. Keep the headline denominator at 93 attempted, not 90 or 73 selected successes.

### Separate local original-MPS diagnostic - v23 passers

All **73 of 73 v23 passing cases** were independently checked locally using the same v23 engine code. `benchmarks/orig_check.py` uses its own MPS parser, maps the returned solution back to original variables, checks row activities and bounds, and compares the recomputed original objective with HiGHS reading the **original MPS file**. This was **not run on Kaggle**, is **not a full-verification gate**, and does **not cover the 20 non-passing cases**.

The per-case pass rule requires HiGHS optimal status, objective error no greater than `max(1e-6, 1e-7 * abs(HiGHS objective))`, relative original-row violation no greater than `1e-6`, and absolute bound violation no greater than `1e-6`. Relative row violation is normalized by `1 + abs(original row RHS)`. The checker applies a 60-second solver limit per case. Observed maxima across the 73 cases: **1.05e-9 relative original-row violation** and **9.3e-10 absolute bound violation**.

**E226 objective convention:** the MPS objective-row RHS entry is **-7.113**, which contributes an additive objective offset of **+7.113**. The solver-reported objective omits this offset: approximately **-18.7519** versus the checker's and original-MPS HiGHS objective **-11.6389**. The corpus/published-reference convention excludes that offset, so E226 remains a pass under that convention; its solution vector is optimal and feasible in the independent check. This disclosure is not a claim that the solver's reported objective includes every original-MPS constant.

To run the checker from the repository root with an existing MPS corpus and `highspy` installed:

```bash
python benchmarks/orig_check.py --engine engine/revised_simplex.py --corpus /path/to/mps_files --out /tmp/orig_check_out.json afiro e226
```

### What changed from the baseline

All **68 v19 passing cases remained passing**. The five gains are BNL1, BNL2, MAROS-R7, MODSZK1 and WOOD1P. Relative to v21, MAROS-R7 is newly passing. The **42/42 regression**, **30/30 synthetic MILP** and **30/30 synthetic QP** checks remained green, with QP maximum objective error **3.197e-13**. Separate selected gates are **EXTRA_NETLIB 7/7** and **CPU_LP 6/7**; CPU_LP is not an all-pass result. These do not erase the 20 named non-passing full-corpus outcomes.

Historical v21: **72/93 passing, 14 timeout, 5 solver_failed, 2 highs_failed**, median **8.4156x over 72 passing pairs**, maximum passing residual **8.15e-10**. Historical v19: **68/93 passing, 16 timeout, 7 solver_failed, 2 highs_failed**, median **12.478x over 68 pairs**, maximum passing residual **9.31e-10**. Historical v18: **57/93 passing, 35 timeout, 1 solver_failed**, median **29.98x over 57 pairs**, maximum passing residual **9.31e-10**. Different passing sets do not make a controlled speedup series.

The full-corpus input source is [the pinned Netlib MPS collection](https://github.com/ozy4dm/lp-data-netlib/tree/56257eea85b433ce6aa67d26156b36385318fd6f/mps_files). The new local and Kaggle ledgers are separate from the older checked-in partial CSVs below. The checked-in engine is synced to the executed v23 solve path, and the parser retains the v19 bulk-construction path. The older checked-in result ledgers below are not the new full-run ledger.

## What is in this repository

- `engine/revised_simplex.py` / `engine/solver_lu_relfeas.py`: matching copies of the executed v23 solve path (product-form basis updates, vectorized pricing, Harris ratio handling, conservative presolve and v19-rule fallback), using NumPy arrays and SciPy LU factorization for linear algebra. SciPy's optimization solvers are **not** called by the core.
- `engine/solver_dantzig_relative.py`: earlier revised-simplex variant retained for provenance, no longer imported by the refinery example and not the engine of record.
- `parsers/mps_fixed.py`: fixed-field MPS parser for rows, RHS, ranges and bounds; it transforms finite lower and upper bounds into nonnegative standard-form variables and extra constraints.
- `parsers/mps_free.py`: extension that splits free variables into positive and negative nonnegative columns. It imports `mps_fixed` from the same directory.
- `benchmarks/run_netlib.py`: example command-line runner to parse an MPS file and report objective, residual, iterations and elapsed time. You supply the input MPS files; they are not included here.
- `results/consolidated_partial_sweep.csv`: older consolidated **partial** sweep ledger, 44 rows (39 passes and five failures). It is not a complete failure inventory. The separate SCSD8 addendum is also counted; the canonical CSV has not yet been regenerated with that row. `results/netlib_all_runs.jsonl` is a prior assembly of the supplied per-case records and is superseded by the CSV.
- `milp/branch_and_bound.py`, `milp/synthetic_3_cases.jsonl`, `milp/property_tests.py`, `milp/property_results.json`: a minimal B&B prototype, three fixed synthetic tests, and 30 randomized small synthetic property tests checked against brute force and SciPy/HiGHS. No MIPLIB coverage.
- `qp/projected_gradient.py`, `qp/synthetic_3_cases.jsonl`, `qp/property_tests.py`, `qp/property_results.json`: positive-definite convex box-bounded projected-gradient QP prototype, three fixed synthetic checks and 30 randomized SPD box-QP synthetic property tests against SciPy. No general linear constraints or public QP benchmark coverage.
- `examples/illustrative_refinery.py`, `examples/illustrative_refinery_result.json`: synthetic refinery LP and retained result, now run on the same LU engine of record as the LP sweep.
- `results/highs_same_machine_14_cases.csv`: a **separate, 14-case** same-machine process-CPU comparison with three measurements per solver and their medians. Do not conflate its denominator with the 40-case LP ledger.

## Older checked-in measurements (historical)

The consolidated partial sweep ledger `results/consolidated_partial_sweep.csv` records 39 distinct passing LP instances and five failed rows. A separate verified SCSD8 addendum (`results/scsd8_verified_addendum.jsonl`) raises that historical distinct verified count to **40**. These are older retained measurements, not the latest local or Kaggle totals above. SHIP04S hit a run limit and is not counted. The retained results match HiGHS objectives and published Netlib reference objectives, with reported maximum primal residual 3.79e-10 across those passes. The five failed rows are not the full inventory of attempts. This is the retained result ledger, not a fresh independent rerun from the checked-in files. The original MPS inputs, reference source links, environment capture and complete run script were not included in the handoff, so the numerical results cannot yet be fully reproduced from this repository alone. Add those artifacts before claiming end-to-end reproducibility.

FIT1D is a separate result file with objective -9146.378092420928 versus HiGHS -9146.378092420926 and published reference -9146.3780924; its own 42.45 CPU seconds versus HiGHS 0.0223 seconds underlines the performance gap.

In that older sweep, attempted cases not counted as verified included SCORPION (singular basis), FINNIS and SCRS8 (parsed HiGHS result mismatches published references), SCAGR25 and MODSZK1 (lost primal feasibility), and reported outside these attached ledgers: DEGEN2, TUFF, SCTAP2, WOOD1P. These historical notes do not override the named outcomes in the newer ledgers above.

The separate 14-case same-machine CSV benchmarked the **earlier pre-LU engine**, not the current LU engine. It reports HiGHS faster on **all 14 tested cases**; the median ratio of our CPU time to HiGHS CPU time across the 14 rows is about **52.7x**. This is an experimental correctness-focused prototype with a substantial performance gap. Timings are process CPU measurements from the supplied CSV, not a universal speed claim. HiGHS is a reference comparator, not part of our solver core.

There is **no GPU LP solve** or general QP solver here. A separate small, pure-integer/binary MILP branch-and-bound prototype is included in `milp/`; its three fixed **synthetic** cases match HiGHS. A seeded property test adds 30/30 randomized small synthetic cases matching brute force and SciPy/HiGHS. These are not MIPLIB coverage or a large-instance MILP guarantee. The separate QP prototype covers only positive-definite convex box-bounded problems; its three fixed synthetic cases and 30/30 seeded randomized SPD box-QP property tests match SciPy objectives. These do not establish general QP or public benchmark coverage. Matrix-operation GPU experiments, if added later, must not be presented as LP solves. The CPU implementation uses sparse basis LU with product-form updates, plus dense fallback. Large refinery models remain out of present verified scope. There is no refinery field data in the checked-in files; any illustrative refinery case should be labelled as synthetic.

## Try the solver on an MPS LP

Install Python 3, NumPy and SciPy, then supply a compatible MPS file:

```bash
python -m pip install numpy scipy
python benchmarks/run_netlib.py path/to/model.mps
```

The runner reports its own result; it does **not** assert a Netlib pass without a separately supplied reference and independently checked model semantics. The parser covers a subset of MPS conventions, and unsupported or malformed cases may fail. Fixed-field formatting and variable bounds matter. Keep the original models and solver environment with any benchmark publication.

## Synthetic refinery example

`examples/illustrative_refinery.py` is a synthetic six-variable refinery LP with two crude streams, capacity, yield, blend-quality and demand constraints. The retained `examples/illustrative_refinery_result.json` reports the same gross margin for the prototype and HiGHS (about $1.19 million/day in its assumed units). Its coefficients are illustrative assumptions, **not refinery operational data**; the figure is not a measured refinery margin or savings estimate.

Run it with `python -m examples.illustrative_refinery` from the repository root after installing NumPy and SciPy.

The repository remains private until its owner chooses otherwise.
