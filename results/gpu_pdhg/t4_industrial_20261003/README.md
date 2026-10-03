# T4 industrial LP measurement

Single-run CUDA PDHG measurements, fp64, shared Colab free Tesla T4 (15 GB), compute capability 7.5; driver 580.82.07, nvcc 13.0.88, two-vCPU host. Engine cap 300s. Reported source: main `58f77af` CPU snapshot plus the then-current `gpu/pdhg.cu`; exact GPU source hash was not included with this CSV.

`total_s` includes parsing. `near_optimal` is approximate, not a certificate of an exact optimum. Recheck uses highspy as a reader on the original expanded MPS. HiGHS objective references are absent for both rail cases. A tolerance in the filename describes the solver stopping rule, not a bound on objective error or every original-model row residual.

At 1e-4, PDS-100 objective relative error is 7.72e-4; FOME21 and PDS-100 original-model relative row violations are 1.23 and 1.02. At 1e-6 they are still 0.0108 and 0.0309, not Netlib-strict feasibility. Both rail cases reach the 300s cap at 1e-6. These measurements do not establish a speed advantage over exact simplex.

[Byte-exact received CSV](ledger.csv).
