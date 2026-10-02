# Round-3 GPU case list (cloud GPU: T4 / A10)

Built from the round-2 evidence on the RTX 5050 laptop in `results/gpu_pdhg/` (`NOTES.md`, `synth_*_tol*.csv`, `CROSSOVER.md`).
`gpu/run_round3.py` reads the table below, so this file is the single source of truth: edit a row here and the harness runs it.
Every row has a round-2 number that justifies it. Round-2 numbers are totals (setup + solve) in seconds on the 5050 laptop;
"cpu10" is the 10-thread CPU PDHG of the same binary. Estimates for the T4 are **not measurements**: 2 x the round-2 GPU total
(a safety factor; the Kaggle T4 solve time for transport at 1M columns was about 1.4 x the local 1e-6 solve time), CPU rows 2 x round-2 cpu10
(cloud vCPUs are fewer and slower). Seeds are the generator defaults, so the models are the round-2 models; where a round-2
HiGHS objective exists the harness fills it in, and a large error against it flags a model mismatch.

Columns: `configs` is `gpu` or `gpu+cpu` (cpu = CPU PDHG with all vCPUs, recorded in the ledger as `threads`);
`check` = independent vector recheck on the original MPS (`gpu/vector_check.py`; `no` where the Python checker would need more RAM than a cloud box has,
and the harness also skips it for 2M+ columns when RAM < 24 GB); `est_s` = estimated total seconds for the row on a T4, all repeats, setup and vector checks included (the Python recheck costs about 31 s per run at 1M columns and 65 s at 2M on the laptop, so it dominates the checked big rows; t100/t316/p10k_e4 rows are measured wall times of the harness on the laptop, the rest are estimates).
Exact CPU simplex is **not** run on the cloud box (no CPU simplex build); the round-2 exact numbers are quoted in the justification where relevant.

## A. Widen the packing margin (GPU was only 3% ahead at 1M columns, 1e-4)
| id | family | size | tol | configs | cap_s | repeats | check | est_s | justification |
|---|---|---|---|---|---|---|---|---|---|
| p1m_e4 | packing | 1000000 | 1e-4 | gpu+cpu | 300 | 3 | yes | 360 | Round 2: GPU 8.99 s vs cpu10 9.31 s (margin 3%, run ranges 8.77-9.08 vs 9.26-9.37); repeat on a different GPU to see if it holds |
| p2m_e4 | packing | 2000000 | 1e-4 | gpu+cpu | 300 | 3 | yes | 635 | Round 2: GPU 18.1 s vs cpu10 22.4 s (cpu1 50.2 s); does the margin widen from 1M to 2M |
| p5m_e4 | packing | 5000000 | 1e-4 | gpu+cpu | 400 | 2 | no | 400 | Extends the 1M -> 2M trend (GPU time 3% below the CPU at 1M, 19% below at 2M); no round-2 number at this size, extrapolated, device memory scales from 0.54 GB at 2M to about 1.4 GB, so it fits a T4; host RAM for parsing the 1.2 GB MPS is the open risk |
| p1m_e6 | packing | 1000000 | 1e-6 | gpu+cpu | 300 | 2 | yes | 454 | Round 2: GPU 30.3 s vs cpu10 52.1 s (GPU 1.7x ahead); confirm on cloud hardware |
| p2m_e6 | packing | 2000000 | 1e-6 | gpu | 300 | 2 | yes | 370 | Round 2: GPU 60.4 s while cpu10 and cpu1 both hit the 120 s limit, and no exact method finishes (cap 300 s here lets the GPU finish with margin on a slower card); CPU baseline not repeated |

## B. Transport where exact CPU gets slow (exact dual simplex 19.5 s at 10M nonzeros in round 2; GPU PDHG was slower there)
| id | family | size | tol | configs | cap_s | repeats | check | est_s | justification |
|---|---|---|---|---|---|---|---|---|---|
| t2000x2500_e6 | transport | 2000x2500 | 1e-6 | gpu | 300 | 1 | no | 300 | Unconverged round-2 cell: GPU hit the 120 s limit at error 2.2e-5 (43,328 iterations, still improving); give it a 300 s cap. Round 2 exact dual: 19.4 s |
| t3000x4000_e4 | transport | 3000x4000 | 1e-4 | gpu | 300 | 1 | no | 300 | 24M nonzeros, 2.4x the largest round-2 transport (GPU 81.1 s at 10M nonzeros, exact dual 19.5 s, CPU PDHG over the limit): device memory about 1.9 GB (scaled from 0.77 GB at 10M nonzeros); finds where the exact dual simplex stops being cheap. Extrapolated GPU time, no round-2 number at this size |

## C. Controls: sizes where the CPU won or the margin was thin, to locate the crossover (gpu+cpu, 3 repeats)
| id | family | size | tol | configs | cap_s | repeats | check | est_s | justification |
|---|---|---|---|---|---|---|---|---|---|
| p10k_e4 | packing | 10000 | 1e-4 | gpu+cpu | 120 | 3 | yes | 8 | Round 2: cpu10 0.050 s beat GPU 0.179 s |
| p100k_e4 | packing | 100000 | 1e-4 | gpu+cpu | 120 | 3 | yes | 33 | Round 2: cpu10 0.593 s beat GPU 0.743 s |
| p300k_e4 | packing | 300000 | 1e-4 | gpu+cpu | 120 | 3 | yes | 90 | Brackets the packing crossover between 100k (CPU ahead) and 1M (GPU 3% ahead) at 1e-4; no round-2 row at this size |
| p10k_e6 | packing | 10000 | 1e-6 | gpu+cpu | 120 | 3 | yes | 10 | Round 2: cpu10 0.226 s beat GPU 0.384 s |
| p100k_e6 | packing | 100000 | 1e-6 | gpu+cpu | 120 | 3 | yes | 58 | Round 2: GPU 2.22 s vs cpu10 2.55 s, a 13% margin, the closest GPU win at 1e-6 |
| t100_e4 | transport | 100x100 | 1e-4 | gpu+cpu | 120 | 3 | yes | 10 | Round 2: cpu10 0.074 s beat GPU 0.187 s |
| t200_e4 | transport | 200x200 | 1e-4 | gpu+cpu | 120 | 3 | yes | 25 | Brackets the transport crossover between 100x100 (CPU ahead) and 316x316 (GPU 1.22 s vs cpu10 3.76 s); no round-2 row at this size |
| t316_e4 | transport | 316x316 | 1e-4 | gpu+cpu | 120 | 3 | yes | 76 | Round 2: GPU 1.22 s vs cpu10 3.76 s; anchors the GPU side of the transport crossover |

Planned total: about 3,130 s at the estimates above (vector checks included) plus model generation (largest: packing 5M, transport 3000x4000, about a minute each on the laptop);
the harness stops starting rows once its budget (default 6,600 s) would be exceeded and records the skipped rows.
Not planned, deliberately: packing 5M at 1e-6 (no round-2 evidence it would finish within a cap that fits the budget), any size at 1e-8, and the Netlib set (round 2 already covers it).
