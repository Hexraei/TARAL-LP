# Mixed-integer benchmark measurement - October 12, 2026

All 65 MIPLIB 3 cases were run in three paired repeats, with a 300-second solver cap per case. Two configurations used the same experimental binary: mixed-integer rounding cuts in both, objective-integrality bound strengthening added to the second. These measurements are not a claim about the unmodified `main` build.

| Reference / counting rule | Control | Objective-integrality strengthening |
| --- | ---: | ---: |
| Catalogue-corrected expected-value agreement | 38/65 | 39/65 |
| Frozen as-run reference agreement | 37/65 | 38/65 |
| Expected-value agreement, finish within 270 seconds in every repeat | 37/65 | 38/65 |

Each row is identical across the three repeats. The only newly solver-optimal case is `vpm1`: control reached the time cap at objective 20 and bound 19.55; strengthened search returned objective/bound 20 in 34.799, 34.821 and 34.767 seconds. `p0548` finishes in both configurations but takes more than 270 seconds in all repeats. `qiu` reaches the reference value without closing the proof and never counts as solved.

## Independent checks and limits

- 390/390 primary slots accounted; 366 returned points independently rechecked on original models. All pass the numerical row/bound/integrality gate of absolute 1e-6 and objective agreement within 1e-6*(1+|objective|). The other 24 slots are capped with no point.
- Every primary cut-verifier result independently replayed. Per configuration: 87 PASS, 69 no-cuts, 27 numerical-INCONCLUSIVE and 12 positive-cuts/no-point INCONCLUSIVE. No strict cut violation was found. Nine cases including `vpm1`, `fiber` and `gesa2` retain numerical-INCONCLUSIVE cut checks; this is not an exact optimality certificate.
- `optimal` plus a passing point and reference agreement is the numerical counting rule. Frozen HiGHS values are not independently proved exact. `cap6000` uses the catalogue value -2451377 in the expected view, versus the historical default-gap value -2451279 in the frozen view. Ten reference-unproved cases have no expected-value stop check.
- Sequential control-first order, two simultaneous single-threaded shards per configuration, three repeats. Fixed order may confound timing. No general speedup claim. `harp2` has a worse capped incumbent in the strengthened configuration.
- Reduced on-host prelaunch tests and a substituted eight-model equivalence set were accepted as a recorded deviation. The planned on-host mutation rerun was not done. Larger development tests and the final independent replays do not retroactively make that prelaunch record complete.

## Measurement identity

Ubuntu 22.04-family host, four-vCPU Xeon Platinum 8272CL, 16 GB RAM; g++ 11.4.0, `-O2 -std=c++17 -pthread`; HiGHS 1.15.1 in the independent audit only.

- Source archive SHA256: `0c32c7a358525718454a1699833262553bbb6c48bb0a2369fa8a343ee9f5548b`
- Solver SHA256: `9774009c472a92c4ce7594929f00a9e08aadb7b1e075fa37216ee1c0d2673d09`
- Evidence archive SHA256: `c03a85bdadf59b0f67f4a7552c2855ac4a71577cbbd03dc726e263619d122417`

The source archive was independently rebuilt twice, reproducing the solver, verifier and source-manifest hashes. This folder contains a judge-facing extraction of the verified [per-case ledger](ledger.csv), not the full private evidence archive or an independently reproducible branch source package. No solver code is merged by this documentation update.
