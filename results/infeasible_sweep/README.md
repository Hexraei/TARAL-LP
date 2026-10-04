# Infeasible-model sweep: evidence archive notes

Describes the evidence behind the tables in `docs/infeasible_sweep.md`. Recorded provenance for two separate builds; neither result supersedes the other.

## What was run
- Tool: `benchmarks/infeasible_sweep.py`, engine plain solve plus `--explain-infeasible`, each explanation checked by
  `benchmarks/infeasibility_explanation_check.py` (HiGHS 1.15.1 as test oracle only).
- Corpus: 29 Netlib infeasible LP models supplied as a package, archive SHA256
  `9e8e1870187a613591a1e41addd2ba6f44162bf9da93402c943017ac8afdf7db`; each file matched the package manifest's `mps_sha256`.
- Engine commit: `ba039383ddea291fa1ff53a33c3493511707031a`, `--time-limit 60`; HiGHS IIS budget 120 s. Native archive binary hash is recorded inside `engine_binary.sha256`. Build for the table in the doc: native g++ 11.4 -O3 -march=native (as of that run).
- Generated LPs: `--generated 60 --seed 7`.

## What is NOT in this repository
- The Netlib corpus. It is not committed and must not be committed: the `gosh` model file contains the problem-statement digits
  as data, and the repository rule bans that literal in tracked content.
- The raw sweep output archive: `sweep_evidence-972a86ec.tgz`, SHA256 `0d187bd83e42943bc6f42d47717200ee7ba9c618b88d18e848310e195ccca927`, recorded October 3, 2026. The archive records an Intel Xeon 2.60GHz Linux 6.1.158+ x86_64 host, g++ 11.4.0, HiGHS 1.15.1. It includes run4/run5 ledgers and run5 explanation witnesses, not the engine binary or corpus.

## Two runs
- Run 1 (native evidence archive, run5 checker at cdedc5b; engine source ba03938): native g++ 11.4 -O3 -march=native: 23 infeasible / 6 numerical_failure.
- Run 2 (Navin's side, independent, generic build): 22 infeasible / 7 numerical_failure (`qual` is `numerical_failure`); `gran` 492 explanation rows vs 1035; `ceria3d` 817 vs 808.
  Generic run: same engine commit `ba039383ddea291fa1ff53a33c3493511707031a`, g++ 11.4.0 `-std=c++17 -O3 -Wall -Wextra -Wpedantic`, no `-march=native`, Intel Xeon 2.60GHz Linux x86_64, engine/explanation budget 60 s. The reviewer independently measured 18 irreducible and 4 reduced_unproven explanations, with all 22 row proofs accepted; cplex1 relaxation failed, qual was not explained. `gran`: 492 rows, all 492 unproven. `ceria3d`: 817 rows, 3 unproven.
- The runs differ. Results depend on the build and the wall-clock budget. Earlier measurements are kept as recorded; nothing was overwritten.

## Claims
- No claim that the engine verdict is certified for any of the 29 models. No "23 certified". HiGHS IIS sizes are reference values only.
- The `cplex1` and `qual` relaxation checks failed in run 1 (engine 0.3826807 vs HiGHS simplex 0.377664455, and 0.00541016 vs 0.00542051); HiGHS interior point gave 0.39088 on `cplex1`.
- Not run: exact-rational certificate verification, MIPLIB infeasible models, repeated wall times, cross-host runs.
