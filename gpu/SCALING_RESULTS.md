# PDHG: matrix-scaling experiments and the opt-in per-row stop test (`--row-rel`)

RTX 5050 laptop, fp64, HiGHS IPM as reference, 300 s cap on the six large cases, single runs (timings about +-3%).
Cases: netlib kennington ken-18, osa-60 (emps-expanded) and Mittelmann pds-100, fome21, rail2586, rail4284.
"Worst row" = max_i viol_i/(1+|b_i|) on the original model, recomputed in numpy from the written x (not by pdhg).

## Summary (read this first)
- **Default behaviour is unchanged.** `--row-rel` defaults to 0 (the old stopping test only). No row of any default-setting ledger improved or was meant to.
- **Matrix scaling is not the lever.** Every scaling variant tried (below) was neutral or worse; the scaling code is unchanged.
- **pds-100 and fome21 reach a small per-row violation at 1e-4 ONLY with `--row-rel 100`.** Without it they stop with worst row 1.02 and 1.23 (objective error 7.7e-4 and 6.7e-5). With it they cost about 3x the time (table 1).
- **`--row-rel 100` is not safe as a default**: it makes three Netlib cases that converge by default fail to converge in 60 s at 1e-4 (table 2).

## Table 1: large cases, default (`--row-rel 0`) vs `--row-rel 100`
Each cell is default -> `--row-rel 100`. Cases where the extra test never binds (osa-60, rail*, at 1e-4) are the same.

| case | tol | iters | total s | objective error | worst row |
|---|---|---|---|---|---|
| pds-100 | 1e-4 | 9344 -> 41152 | 8.9 -> 30.7 | 7.7e-04 -> 2.9e-07 | 1.0e+00 -> 5.1e-03 |
| fome21 | 1e-4 | 7232 -> 23104 | 2.1 -> 5.7 | 6.7e-05 -> 1.5e-07 | 1.2e+00 -> 1.6e-03 |
| ken-18 | 1e-4 | 14464 -> 16512 | 2.9 -> 3.3 | 5.1e-08 -> 9.4e-09 | 9.7e-03 -> 1.3e-03 |
| osa-60 | 1e-4 | 3072 -> 3072 | 3.9 -> 3.9 | 7.7e-05 -> 7.7e-05 | 6.0e-04 -> 6.0e-04 |
| rail2586 | 1e-4 | 36928 -> 36928 | 99.0 -> 99.5 | 5.1e-06 -> 5.1e-06 | 9.7e-04 -> 9.7e-04 |
| rail4284 | 1e-4 | 34368 -> 34368 | 114.5 -> 114.8 | 6.0e-06 -> 6.0e-06 | 1.2e-03 -> 1.2e-03 |
| pds-100 | 1e-6 | 35840 -> 46336 | 26.1 -> 33.5 | 2.0e-06 -> 4.2e-09 | 3.1e-02 -> 3.6e-05 |
| fome21 | 1e-6 | 15552 -> 34624 | 4.0 -> 8.3 | 1.7e-06 -> 4.6e-10 | 1.1e-02 -> 2.1e-05 |
| ken-18 | 1e-6 | 43200 -> 58304 | 7.8 -> 10.4 | 3.4e-08 -> 9.9e-11 | 5.6e-04 -> 6.2e-06 |
| osa-60 | 1e-6 | 9792 -> 21952 | 10.6 -> 22.5 | 1.6e-06 -> 1.1e-06 | 7.0e-05 -> 4.6e-06 |
| rail2586 | 1e-6 | 113472 -> 113408 | 300.1 -> 300.0 (time limit, both) | 2.1e-07 -> 2.1e-07 | 4.9e-05 -> 4.8e-05 |
| rail4284 | 1e-6 | 91776 -> 91840 | 300.1 -> 300.1 (time limit, both) | 2.2e-08 -> 2.3e-08 | 1.3e-04 -> 1.3e-04 |

Cost of `--row-rel 100`: pds-100 at 1e-4 8.9 s -> 30.7 s and fome21 2.1 s -> 5.7 s (about 3x); at 1e-6 the four non-rail cases take 1.3-2.1x longer for tighter answers. The rail cases hit the cap either way.

## Table 2: regressions at `--row-rel 100` (Netlib corpus, 1e-4, 60 s limit; committed ledger = default behaviour)
| case | default (= committed ledger) | `--row-rel 100` |
|---|---|---|
| bnl1 | near_optimal, 6,083,008 iters, 36.9 s | time limit, 9,892,032 iters; objective error 5.7e-4 -> 1.2e-3 |
| perold | near_optimal, 4,039,360 iters, 27.1 s | time limit, 8,785,280 iters |
| pilot.we | near_optimal, 201,472 iters, 1.7 s | time limit, 7,514,176 iters, not converged |

Cause: the extra criterion needs the per-row-relative residual <= 1e-2 (K x tol); these models sit at 5e-4..7e-4 row residual while the gap is already small, and the iteration never gets below. With `--row-rel 0` all three reproduce the committed iteration counts exactly. The other 91 of 94 files had no status change at `--row-rel 100`, and the recheck found no vector failures (worst recheck row violation 0.038, was 130 in the committed ledger).

## Why the stop test
The old test is ||viol||_2 <= tol (1 + ||b||_2) over the original rows, so with a large ||b|| one row can be off by about its own size and the test still passes. The added test is ||viol_i / (1 + |b_i|)||_2 <= K tol (K = `--row-rel`, 0 = off). K tried: 1 (pds-100 3.8x slower, osa-60 objective error 1.6e-4), 30 (osa-60 objective error 3.6e-4), 100 (adopted as the opt-in value).

## Scaling variants tried (1e-4 on pds-100, fome21, ken-18, osa-60; baseline Ruiz 10 passes inf-norm + Pock-Chambolle)
| variant | verdict |
|---|---|
| Ruiz 20 / 50 passes, or 0 passes | identical to baseline on pds/fome/ken (already converged); osa-60 +-1 digit. No gain. |
| no Pock-Chambolle | pds-100 4x slower, ken-18 8x slower, osa-60 objective error 1e-1. Worse. |
| L2 Ruiz (10 / 20 passes) | pds objective error 4.7e-5 / 6.9e-6 but worst row 2.3 / 2.7, 1.5-3x slower, osa-60 slower. Worse. |
| rhs-aware Ruiz (|b| in the row measure, weight 1 / 0.1) | pds worst row 0.058 / 0.055 but 14x / 3.5x slower, ken-18 300 s timeout / 114 s, fome21 objective error 2e-4. Worse. |
| rhs-aware, no Pock-Chambolle | every case hits the 300 s cap. Worse. |

## Determinism notes
- GPU, default binary vs the previous binary: `.sol` and `.dual` byte-identical on fome21 and ken-18 (1e-4).
- CPU: for pilot.we (1 thread, 1e-4) the iteration count is identical (201,472) but the written solution is not byte-identical to the previous binary. The baseline binary reproduces itself, so the cause is the new binary: floating-point reordering from the extra reduction term (computed even when `--row-rel` is 0). Not a behavioural change at the stopping level.

## Validation of the default setting (no flags), Oct 3
- Full 94-file sweeps at 1e-4 and 1e-6 (`gpu/bench_pdhg.py`, 3 repeats, `--repeat-all`, 60 s limit, vector recheck of every run on the original MPS) compared row by row with `results/gpu_pdhg/netlib_94files_tol1e-{4,6}.csv`: no status change, no objective-error regression, no recheck failure. Iterations differ only on the wall-clock time-limit rows (2 at 1e-4: greenbea, pilot.ja; 8 at 1e-6: bnl1, bnl2, greenbea, greenbeb, perold, pilot, pilot.ja, pilot.we), which differ across repeats by construction. Worst recheck row violation unchanged (130 at 1e-4, 2.4 at 1e-6, as in the committed ledgers).
- CPU, 1 thread, 1e-4: iteration counts identical to the previous binary on adlittle, afiro, kb2 and pilot.we; solutions byte-identical except pilot.we (see above).
