# Maros-Meszaros ledger: reading the 5 objective mismatches (local re-run, current build)

`ledger.csv` is unchanged (cloud container, pass rule unchanged: objective within 1e-6 of the HiGHS reference on the same file).
Five rows fail that rule although the solver reported optimal. A local re-run of those five reproduces the same numbers; each was then
checked against something other than the default HiGHS answer:

| Instance | taral | HiGHS default reference | Finding |
| --- | --- | --- | --- |
| DPKLO1 | 0.370096217 | 0.712522210 | HiGHS reads the RHS section of this file wrongly (the RHS set name is the numeric string `1`, so it drops almost all right-hand sides: 1 nonzero row bound instead of 62). Re-writing only the RHS set name to `RHS` makes HiGHS return 0.370096217114, equal to taral. The independent checker (own parser) finds row violation 2e-14 for taral's point. The 36.5 "row violation" in the ledger comes from checking against HiGHS's misread model. |
| QBORE3D | 3100.20080 | 3102.13880 | HiGHS with `qp_regularization_value=1e-12` returns 3100.20080176; taral's dual bound 3100.20080146 (gap 1.4e-10). The default HiGHS answer is the less accurate one. |
| GOULDQP2 | 1.84277621e-4 | 1.88202517e-4 | HiGHS with regularization 1e-12 returns 1.8427450e-4; taral's dual bound 1.8427056e-4 (gap 7e-9). |
| PRIMALC8 | -18309.42979 | -18309.26684 | HiGHS with regularization 1e-12 returns -18309.4297884; taral's dual bound -18309.4297895. |
| QFFFFF80 | 873147.46057 | 873149.27707 | taral is lower, feasible (violation 3e-11) and self-reported gap 1.3e-10. HiGHS with regularization 1e-12 returned 0.0 (unusable), so no second reference exists. Not independently certified. |

**Summary: on DPKLO1, QBORE3D, GOULDQP2 and PRIMALC8 the taral answer is the accurate side (HiGHS's default reference is wrong or less accurate); QFFFFF80 stays uncertified.**

Caveats: the dual bounds above are taral's own numbers (not recomputed by independent code); the rows stay FAIL under the stated rule.
No solver change was needed for these five. Not run here: the 4 numerical failures. The 18 rows without a HiGHS reference are covered in the next section.

## The 18 rows without a HiGHS reference, rechecked with `qp_regularization_value=1e-12`

`benchmarks/qp_recheck_unreferenced.py` (HiGHS 1.15.1, 300 s per instance, same files) writes `recheck_unreferenced.csv`.
Outcome: **2 of 18 get a reference and both agree** with TARAL under the ledger rule (relative error <= 1e-6, violation <= 1e-6):
QCAPRI (error 1.1e-12; the default-HiGHS reference had failed with a solve error) and QSCTAP3 (error 8.5e-10).
**The other 16 stay unreferenced**: HiGHS again returns no optimal answer (not set, solve error, time limit, or, for QSTAIR,
"unbounded" while TARAL reports a feasible optimum with violation 7.9e-9). For these 16 the TARAL answers rest on their own
recheck only (independent row/bound violation and objective recomputation in the ledger; all relative violations <= 8e-9), not on
agreement with a second solver. The ledger rows and verdicts are unchanged (the pass rule needs a reference).
HiGHS's status on a given file is not stable across settings and runs (QCAPRI: solve error by default, optimal here), so these
rows are re-checkable but should not be read as evidence about TARAL either way.
