# Second-host reproduction and 60-second read - October 16, 2026

Backing records for the two additions made to the repository README on October 16, 2026: the second-host reproduction of the 42/65 expected-value count, and the 60-second read (28/65 vs 23/65). Counting rules and reference views are those of the [October 13 measurement](../milp_paired_20261013/).

## Second-host reproduction, 300 seconds per case (`ledger_r16b.csv`)

All 65 MIPLIB 3 cases ran in three paired repeats at 300 seconds per case on a second host (Azure D4s_v4, 4 vCPUs, Intel Xeon Platinum 8272CL, Ubuntu 22.04) on October 15, 2026. Two configurations on the same binary: defaults (control) and the measured opt-in configuration (`TARAL_MIR=1 TARAL_OBJINT=1 TARAL_MIR_AGG=1 TARAL_DIVE_ROUNDS=1`; single-threaded BLAS). Binary sha256 `6a6ac74fb30c0fa09554cbe8ed6671017b6d31d95236ed0981d39164470d11f4` (deterministic double build from source archive `f871c95efe966e15c69500846f86e523dea83c976b0150a9810e149a92900558`, build script `f681efb8e68a7fe8a3b7f8ae682f72e0de30807ed450465466c1731c4d5820cf`; g++ 11.4.0, `-O2 -std=c++17 -pthread`).

| Counting rule | Control (defaults) | Measured opt-in |
|---|---|---|
| Frozen reference | 33/65 | 41/65 |
| Expected value | 34/65 | 42/65 |
| Robust expected value (`p0548` near-cap exclusion) | 34/65 | 41/65 |

Counts are identical in all three repeats. The measured opt-in configuration converts `fiber`, `gesa2`, `gesa3`, `pp08a`, `gesa2_o`, `p0548`, `pp08aCUTS` and `vpm1` over defaults, with no reversals; this reproduces the October 13 expected-value count of 42/65 on the second host. 390/390 case runs valid. Every returned point passes the independent audit on the original models (62/62 audited runs per arm-repeat block; tolerances 1e-6). Per-case rows: `ledger_r16b.csv`. Time-limited counts depend on host speed; default-configuration counts on this host are not comparable to the historical host.

## Corroborating run, second host, 300 seconds per case (`ledger_r16c.csv`)

A separate paired run on another Azure D4s_v4 host (pair300cvm: 4 vCPUs, Intel Xeon Platinum 8272CL @ 2.60GHz, Ubuntu 22.04, kernel 6.8.0-1070-azure) on October 16, 2026 reproduced the same counts with a different build. All 65 MIPLIB 3 cases ran in three repeats at 300 seconds per case under the same measured opt-in configuration (`TARAL_MIR=1 TARAL_OBJINT=1 TARAL_MIR_AGG=1 TARAL_DIVE_ROUNDS=1`; single-threaded BLAS; CPU affinity 0-3); all other experimental features were off.

| Counting rule | Repeat 1 | Repeat 2 | Repeat 3 |
|---|---|---|---|
| Frozen reference | 41/65 | 41/65 | 41/65 |
| Expected value | 42/65 | 42/65 | 42/65 |
| Robust expected value (`p0548` excluded) | 41/65 | 41/65 | 41/65 |

195/195 case runs valid. Every returned point passes the independent audit on the original models (all tolerances 1e-6; maximum absolute violation 0.0). Per-case rows: `ledger_r16c.csv` (same column semantics as the October 13 ledger, with the robust expected-value rule as an added column). Binary sha256 `a676a4ff8a7b05add0dda06695ee867ec7a21c7e3988433bd4b32151dfab37dd` (g++ 11.4.0, `-O2 -std=c++17 -pthread`).

## 60-second read (`ledger_60s_cutoff_read.csv`)

Reading the October 15 second-host runs above at a 60-second cutoff: a case counts when the solver reached `optimal` status within 60 seconds of wall time. This is a cutoff read of the 300-second runs, not a separate 60-second-capped run.

| Cutoff rule | Control (defaults) | Measured opt-in |
|---|---|---|
| Proven within 60 s | 23/65 | 28/65 |
| Proven within 54 s (robust rule) | 23/65 | 26-27/65 |

The 60-second counts are identical in all three repeats, with identical case sets: the measured opt-in configuration adds `fiber`, `gesa3`, `gesa3_o`, `pp08a`, `pp08aCUTS` and `vpm1` over defaults, and defaults add `dcmulti`. The robust-rule count is 27/65 in repeats 2 and 3 and 26/65 in repeat 1 (`blend2` proven in 54.2 s). Per-case rows: `ledger_60s_cutoff_read.csv` (derived from `ledger_r16b.csv`; same binary, host and configuration). Cutoff counts depend on host speed; the 300-second counts above are the stable reference.
