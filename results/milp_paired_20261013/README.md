# Mixed-integer benchmark measurement - October 13, 2026

All 65 MIPLIB 3 cases were run in three paired repeats, with a 300-second solver cap per case. Two configurations used the same experimental binary: mixed-integer rounding cuts and objective-integrality bound strengthening in both, multi-row aggregation added to the second. This is an experimental snapshot, not the unmodified `main` build. The [October 12 measurement](../milp_paired_20261012/) is historical.

| Reference / counting rule | Control | Row aggregation |
| --- | ---: | ---: |
| Catalogue-corrected expected-value agreement | 39/65 | 42/65 |
| Frozen as-run reference agreement | 38/65 | 41/65 |
| Expected-value agreement, finish within 270 seconds in every repeat | 38/65 | 41/65 |

Each row is identical across the three repeats. Only `pp08a`, `gesa2_o` and `pp08aCUTS` convert from capped to solver-optimal status; no solved case reverses.

| Case | Control outcome | Aggregation objective | Aggregation wall seconds, three repeats |
| --- | --- | ---: | --- |
| pp08a | Capped, incumbent about 7650 | about 7350 | 3.787 / 3.780 / 3.775 |
| gesa2_o | Capped, no returned point | 25779856.371697895 | 94.292 / 94.326 / 93.847 |
| pp08aCUTS | Capped, incumbent 7760 | 7350 | 3.340 / 3.358 / 3.385 |

The conversions use the solver's existing relative-gap semantics. `gesa2_o` has final bound 25779832.924077682, below its objective; `pp08aCUTS` also has a small numerical gap. All three retain numerical-INCONCLUSIVE cut checks. `p0548` solves in both configurations at 282-289 seconds and is the sole solved case excluded by the 270-second margin. `qiu` reaches its reference value without closing the search and is never counted as solved.

## Independent checks and limits

- 390/390 primary slots valid, zero missing or invalid slots and no expected-value defects. All 369 returned points independently pass original-model row/bound/integrality checks at absolute 1e-6 and objective agreement within 1e-6*(1+|objective|). The other 21 are capped with no point.
- All 390 primary cut-verifier RESULT lines independently replayed with exact archived agreement and no strict violations across 30,279 emitted cuts. Control: 87 PASS, 69 no-cuts, 27 numerical-INCONCLUSIVE, 12 positive-cuts/no-point INCONCLUSIVE. Aggregation: 69, 63, 54, 9 respectively. The 24 calibration checks independently replay no-cuts. Numerical-INCONCLUSIVE remains a separate axis; absence of a strict violation is not exact validation of all cuts at an exact source-feasible vector.
- The numerical count requires solver `optimal` status, a passing point audit and reference agreement. It is not exact proof. `cap6000` uses catalogue value -2451377 in the expected view versus historical default-gap value -2451279 in the frozen view. Ten reference-unproved cases have no expected-value stop check.
- Sequential control-first order, two simultaneous single-threaded shards per configuration, three repeats. Fixed order may confound timing. No broad speedup claim: median wall-time ratio aggregation/control is about 1.0041 across the 39 cases solved in all six slots. `vpm1` slows from about 34.6 to 56 seconds. The three conversions remain useful measured outcomes, not evidence of a general speed advantage.
- No block or within-shard overlap; 28,466 process samples show at most two tracked solvers and no outside solver process. Environment observation is unverified for 21 very short slots, with no observed mismatch; trace configuration attestations remain available.
- Knob-unset prelaunch comparison: **equivalence 59/65, 6 timing-consistent differences accepted**. This does not prove all-case equivalence. The source package was revised before launch but reached independent review after launch; identity was reconciled during the run, with solver source/build unchanged from the reviewed version. This routing deviation remains recorded.

## Measurement identity

Ubuntu 22.04-family host, four-vCPU Xeon Platinum 8272CL, 16 GB RAM; g++ 11.4.0, `-O2 -std=c++17 -pthread`; HiGHS 1.15.1 in the independent audit only.

- Source archive SHA256: `f871c95efe966e15c69500846f86e523dea83c976b0150a9810e149a92900558`
- Solver SHA256: `6a6ac74fb30c0fa09554cbe8ed6671017b6d31d95236ed0981d39164470d11f4`
- Evidence archive SHA256: `e795b5c02cd87c294b316133e86bf4fbcb11ac684ae1083075156f466a0f9247`

Two independent archive rebuilds reproduce the solver, verifier and source-manifest hashes. This folder contains a judge-facing extraction of the verified [per-case ledger](ledger.csv), not the full private evidence archive or an independently reproducible branch source package. No solver code is merged by this documentation update.
