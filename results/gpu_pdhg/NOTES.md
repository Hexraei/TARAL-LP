# GPU PDHG measurements (local, measured; prototype, near-optimal)

Hardware: RTX 5050 Laptop GPU (compute capability 12.0, 8 GB), Intel Core 7 240H (10 cores / 16 threads), CUDA 12.9.
Solver: `gpu/pdhg.cu` (CUDA runtime and C++17 standard library only; own kernels; matrix resident on the GPU) with a
CPU backend of the same algorithm. PDHG is a near-optimal method: every result is reported with its objective error
against HiGHS on the original MPS file, never as an exact solve, and is a separate line from the simplex headline.

Files
- `netlib.csv/json`: seven Netlib cases, tol 1e-4, median of 3 runs (re-measured Oct 2), GPU versus 1 and 10 CPU threads (+ exact simplex time).
- `synth_transport_tol{1e-4,1e-6}.csv/json`, `synth_packing_tol{1e-4,1e-6}.csv/json`: round 2 (Oct 2, current build, vector checks, 3 repeats of every config, exact primal and dual simplex rows); tables below.
- Round 1, kept for the record (earlier simplex build, no vector checks): `synth_transport.csv/json`, `synth_packing_part1.log`, `packing_2000000_gpu.json`: synthetic instances from `gpu/gen_lp.py`
  (fixed seed; transport and packing LPs; clearly synthetic, not industrial data).
- `netlib_94files_tol1e-4.csv`, `netlib_94files_tol1e-6.csv`: correctness sweep over all 94 corpus files (the 93 cases plus TRUSS, one row each) at two tolerances
  (objective error versus HiGHS per case). At 1e-4: 92 near-optimal, 2 time limits (greenbea, pilot.ja); at 1e-6: 86 and 8.
  Re-measured Oct 2 with the vector-check harness (every finished run rechecked on the original MPS, 3 repeats of every
  config, `--repeat-all`; no recheck failure in any row). Versus the first sweep (88 and 6 at 1e-6) greenbeb and perold
  now hit the 60 s limit: both sit at the edge, so the near-optimal counts at 1e-6 are not stable. Time-limited rows differ
  in iterations across repeats by construction (wall-clock stop). At 1e-4 the worst objective error among the
  near-optimal rows is 4.8e-2 (forplan; nine rows exceed 1e-3): the stopping test is scaled, read the error column.

Measured crossover (GPU total time including setup and transfer versus the best CPU configuration)
- Netlib cases tested (up to 73k nonzeros): the 10-thread CPU is faster on every case.
- Transport family: GPU wins from about 2e5 nonzeros (transport_316x316: 1.27 s versus 3.86 s); at 2 million nonzeros 15.5 s
  versus 76 s; at 10 million nonzeros 82 s while the CPU configurations hit the 120 s limit (errors 2.2e-3 and 4.2e-4).
- Packing family: the 10-thread CPU wins at 100k columns (0.56 s versus 0.65 s); the GPU wins at 1 million (7.35 s versus 9.35 s).
- The crossover depends on the instance, not only on nonzeros. No claim is made beyond these instances.

Caveats
- (Round 1 caveat; the round 2 section below reports the measured dimensions of every instance.) Dimensions of the packing instances other than the 2 million column case were inferred, not measured; packing_1M has no
  HiGHS reference and packing_2M no reference and a single GPU run (HiGHS did not finish; the CPU runs wrote no output).
- Convergence at the stated tolerance can hide a large objective error on hard cases; always read the error column.
- Timings are local, on a laptop that was also running other work (load average 1.2 to 3.5); not Kaggle, no speed claim
  versus HiGHS or commercial solvers. GPU first runs after idle can be slower while the card leaves a low-power state.

## Synthetic families, round 2 (Oct 2, RTX 5050 Laptop + Intel Core 7 240H, current main build 5e93146; supersedes the round-1 synthetic numbers above)
Protocol: `gpu/bench_pdhg.py --configs gpu,cpu1,cpu10 --repeats 3 --repeat-all --time-limit 120 --simplex --simplex-methods simplex,dual --max-load 2` at tol 1e-4 and 1e-6;
every finished run's written primal and dual vectors are rechecked on the original MPS by `gpu/vector_check.py` (0 recheck failures in all 48 PDHG rows; the exact rows are single runs).
"limit (T)" = stopped at the 120 s limit without an answer (T = wall seconds). Exact simplex rows are the C++ engine, one thread; where an exact solve finished its
objective matched HiGHS (worst relative difference 4.8e-15). GPU objective error is against HiGHS; packing_1M and packing_2M have no HiGHS reference (HiGHS did not finish), so no error is shown.
Load average during the runs was 1.0 to 2.6 for single-thread configs (the 10-thread CPU configs raise it themselves).

### transport, PDHG tol 1e-4 (time limit 120 s; median total of 3 repeats, seconds; objective error vs HiGHS where a reference exists)

| Case | rows x cols (nnz) | GPU PDHG | CPU 10 threads | CPU 1 thread | exact primal simplex | exact dual simplex | GPU objective error |
|---|---|---|---|---|---|---|---|
| transport_100x100 | 200 x 10000 (20000) | 0.187 | 0.0739 | 0.188 | 0.267 | 0.0146 | 4.99e-05 |
| transport_316x316 | 632 x 99856 (199712) | 1.22 | 3.76 | 10.1 | 18.7 | 0.17 | 2.67e-05 |
| transport_1000x1000 | 2000 x 1000000 (2000000) | 15.5 | 74.6 | 118 | limit (120) | 2.88 | 9.58e-05 |
| transport_2000x2500 | 4500 x 5000000 (10000000) | 81.1 | limit (121) | limit (121) | limit (120) | 19.5 | 3.49e-05 |

### transport, PDHG tol 1e-6 (time limit 120 s; median total of 3 repeats, seconds; objective error vs HiGHS where a reference exists)

| Case | rows x cols (nnz) | GPU PDHG | CPU 10 threads | CPU 1 thread | exact primal simplex | exact dual simplex | GPU objective error |
|---|---|---|---|---|---|---|---|
| transport_100x100 | 200 x 10000 (20000) | 0.225 | 0.123 | 0.33 | 0.289 | 0.015 | 6.23e-07 |
| transport_316x316 | 632 x 99856 (199712) | 1.94 | 6.37 | 17.8 | 19 | 0.182 | 7.98e-07 |
| transport_1000x1000 | 2000 x 1000000 (2000000) | 37.4 | limit (120) | limit (120) | limit (120) | 2.59 | 2.52e-06 |
| transport_2000x2500 | 4500 x 5000000 (10000000) | limit (120) | limit (121) | limit (121) | limit (120) | 19.4 | 2.18e-05 |

### packing, PDHG tol 1e-4 (time limit 120 s; median total of 3 repeats, seconds; objective error vs HiGHS where a reference exists)

| Case | rows x cols (nnz) | GPU PDHG | CPU 10 threads | CPU 1 thread | exact primal simplex | exact dual simplex | GPU objective error |
|---|---|---|---|---|---|---|---|
| packing_10000 | 5000 x 10000 (49980) | 0.179 | 0.0503 | 0.104 | limit (120) | 68 | 1.14e-06 |
| packing_100000 | 50000 x 100000 (499978) | 0.743 | 0.593 | 1.42 | limit (120) | limit (120) | 1.01e-06 |
| packing_1000000 | 500000 x 1000000 (4999980) | 8.99 | 9.31 | 19.8 | limit (120) | limit (120) | - |
| packing_2000000 | 1000000 x 2000000 (9999978) | 18.1 | 22.4 | 50.2 | limit (120) | limit (120) | - |

### packing, PDHG tol 1e-6 (time limit 120 s; median total of 3 repeats, seconds; objective error vs HiGHS where a reference exists)

| Case | rows x cols (nnz) | GPU PDHG | CPU 10 threads | CPU 1 thread | exact primal simplex | exact dual simplex | GPU objective error |
|---|---|---|---|---|---|---|---|
| packing_10000 | 5000 x 10000 (49980) | 0.384 | 0.226 | 0.797 | limit (120) | 64.4 | 5.39e-09 |
| packing_100000 | 50000 x 100000 (499978) | 2.22 | 2.55 | 9.47 | limit (120) | limit (120) | 1.58e-09 |
| packing_1000000 | 500000 x 1000000 (4999980) | 30.3 | 52.1 | limit (121) | limit (120) | limit (120) | - |
| packing_2000000 | 1000000 x 2000000 (9999978) | 60.4 | limit (121) | limit (123) | limit (120) | limit (120) | - |

What changed against round 1 and what did not
- Transport: the exact dual simplex (new in this round's table) is faster than every PDHG configuration at every size: 2.9 s against 15.5 s (GPU) at 2 million nonzeros, 19.5 s against 81 s at 10 million.
  The exact primal simplex does not finish from 1000 x 1000 up. Round 1 compared the GPU only with the primal simplex; against the dual simplex the GPU PDHG has no advantage on this family.
  GPU against the 10-thread CPU at 1e-4 reproduces round 1 (1.22 s against 3.76 s at 316 x 316; 15.5 s against 74.6 s at 1000 x 1000). At 1e-6, transport_2000x2500 does not converge on the GPU in 120 s (error 2.2e-5).
- Packing: the exact simplex methods finish only on the 10,000-column case (dual 64 to 68 s, primal not in 120 s) and no exact method finishes from 100,000 columns up, so PDHG is the only route to an answer there.
  GPU against the 10-thread CPU: at 1e-4 the CPU is faster at 10k and 100k columns, the GPU only marginally at 1M (8.99 s against 9.31 s; round 1 had 7.35 against 9.35, so that margin shrank and is within what a laptop shows from run to run) and clearly at 2M (18.1 s against 22.4 s). At 1e-6 the GPU wins from 100k columns (2.22 s against 2.55 s; 30.3 s against 52.1 s at 1M; 60.4 s while both CPU configurations hit the limit at 2M), and the CPU wins at 10k.
- Packing_2000000 now has CPU rows that finish at 1e-4 (cpu10 22.4 s, cpu1 50.2 s); round 1 recorded that the CPU runs wrote no output there.
- The simplex on packing_100000 is the same engine whose warm-started variant is described in CROSSOVER.md: still no exact answer in 120 s.
