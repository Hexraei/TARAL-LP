# TARAL-LP

An experimental, inspectable CPU linear-programming prototype for this submission. The current core is a two-phase revised simplex solver for small LPs, with an MPS parser path. It is not a production refinery planner or a replacement for HiGHS, CPLEX, or Xpress.

## What is in this repository

- `engine/revised_simplex.py` / `engine/solver_lu_relfeas.py`: identical copies of the LU-based engine of record for the 25-case regression and wider sweeps, using NumPy arrays and SciPy LU factorization for linear algebra. SciPy's optimization solvers are **not** called by the core.
- `engine/solver_dantzig_relative.py`: earlier revised-simplex variant used for the illustrative refinery case, not the engine of record for the 39-case LP ledger.
- `parsers/mps_fixed.py`: fixed-field MPS parser for rows, RHS, ranges and bounds; it transforms finite lower and upper bounds into nonnegative standard-form variables and extra constraints.
- `parsers/mps_free.py`: extension that splits free variables into positive and negative nonnegative columns. It imports `mps_fixed` from the same directory.
- `benchmarks/run_netlib.py`: example command-line runner to parse an MPS file and report objective, residual, iterations and elapsed time. You supply the input MPS files; they are not included here.
- `results/netlib_all_runs.jsonl`: 44 retained attempt records, with source-ledger filename preserved per row; 39 distinct passes and five non-counted failures. The initial 25-case and wider-sweep ledgers were consolidated into this single checked-in file.
- `milp/branch_and_bound.py`, `milp/synthetic_3_cases.jsonl`: a minimal B&B prototype and three synthetic test records, not MIPLIB results.
- `examples/illustrative_refinery.py`, `examples/illustrative_refinery_result.json`: synthetic refinery LP and retained result.
- `results/highs_same_machine_14_cases.csv`: a **separate, 14-case** same-machine process-CPU comparison with three measurements per solver and their medians. Do not conflate its denominator with the 39-case LP ledger.

## Measured status

The supplied result records, consolidated in `results/netlib_all_runs.jsonl`, show **39 distinct passing LP instances**: 25/25 selected cases in the first ledger and 14 additional passes in wider-sweep and separate-result files. The retained results match HiGHS objectives and published Netlib reference objectives, with reported maximum primal residual 3.79e-10 across those passes. This is the retained result ledger, not a fresh independent rerun from the checked-in files. The original MPS inputs, reference source links, environment capture and complete run script were not included in the handoff, so the numerical results cannot yet be fully reproduced from this repository alone. Add those artifacts before claiming end-to-end reproducibility.

FIT1D is a separate result file with objective -9146.378092420928 versus HiGHS -9146.378092420926 and published reference -9146.3780924; its own 42.45 CPU seconds versus HiGHS 0.0223 seconds underlines the performance gap.

Other attempted cases not counted as verified include SCORPION (singular basis), FINNIS and SCRS8 (parsed HiGHS result mismatches published references), SCAGR25 and MODSZK1 (lost primal feasibility), and reported outside these attached ledgers: DEGEN2, TUFF, SCTAP2, WOOD1P. Those further failures need their own raw records before being independently audited.

The separate 14-case same-machine CSV reports HiGHS faster on **all 14 tested cases**; the median ratio of our CPU time to HiGHS CPU time across the 14 rows is about **52.7x**. This is an experimental correctness-focused prototype with a substantial performance gap. Timings are process CPU measurements from the supplied CSV, not a universal speed claim. HiGHS is a reference comparator, not part of our solver core.

There is **no GPU LP solve** or QP solver here. A separate small, pure-integer/binary MILP branch-and-bound prototype is included in `milp/`; its three **synthetic** cases match HiGHS in the attached ledger. That is not MIPLIB coverage and not a large-instance MILP guarantee. Matrix-operation GPU experiments, if added later, must not be presented as LP solves. The CPU implementation densifies the model and refactorizes a basis, so large refinery models are out of present scope. There is no refinery field data in the checked-in files; any illustrative refinery case should be labelled as synthetic.

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

## Next work

Keep matrices sparse, broaden verified problem coverage, add a fully reproducible harness with inputs and reference provenance, build and measure a GPU LP path, then broaden MILP beyond three synthetic tests and tackle QP separately. These planned items are not counted as implemented.

No license has been selected yet. The repository remains private until its owner chooses otherwise.
