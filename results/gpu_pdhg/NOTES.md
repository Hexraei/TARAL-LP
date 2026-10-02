# GPU PDHG measurements (local, measured; prototype, near-optimal)

Hardware: RTX 5050 Laptop GPU (compute capability 12.0, 8 GB), Intel Core 7 240H (10 cores / 16 threads), CUDA 12.9.
Solver: `gpu/pdhg.cu` (CUDA runtime and C++17 standard library only; own kernels; matrix resident on the GPU) with a
CPU backend of the same algorithm. PDHG is a near-optimal method: every result is reported with its objective error
against HiGHS on the original MPS file, never as an exact solve, and is a separate line from the simplex headline.

Files
- `netlib.csv/json`: seven Netlib cases, tol 1e-4, median of 3 runs, GPU versus 1 and 10 CPU threads (+ exact simplex time).
- `synth_transport.csv/json`, `synth_packing_part1.log`, `packing_2000000_gpu.json`: synthetic instances from `gpu/gen_lp.py`
  (fixed seed; transport and packing LPs; clearly synthetic, not industrial data).
- `netlib_all93_tol1e-4.csv`, `netlib_all93_tol1e-6.csv`: correctness sweep over all Netlib cases at two tolerances
  (objective error versus HiGHS per case). At 1e-4: 92 near-optimal, 2 time limits; at 1e-6: 88 and 6.

Measured crossover (GPU total time including setup and transfer versus the best CPU configuration)
- Netlib cases tested (up to 73k nonzeros): the 10-thread CPU is faster on every case.
- Transport family: GPU wins from about 2e5 nonzeros (transport_316x316: 1.27 s versus 3.86 s); at 2 million nonzeros 15.5 s
  versus 76 s; at 10 million nonzeros 82 s while the CPU configurations hit the 120 s limit (errors 2.2e-3 and 4.2e-4).
- Packing family: the 10-thread CPU wins at 100k columns (0.56 s versus 0.65 s); the GPU wins at 1 million (7.35 s versus 9.35 s).
- The crossover depends on the instance, not only on nonzeros. No claim is made beyond these instances.

Caveats
- Dimensions of the packing instances other than the 2 million column case are inferred, not measured; packing_1M has no
  HiGHS reference and packing_2M no reference and a single GPU run (HiGHS did not finish; the CPU runs wrote no output).
- Convergence at the stated tolerance can hide a large objective error on hard cases; always read the error column.
- Timings are local, on a laptop that was also running other work (load average 1.2 to 3.5); not Kaggle, no speed claim
  versus HiGHS or commercial solvers. GPU first runs after idle can be slower while the card leaves a low-power state.
