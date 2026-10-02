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

Caveats: the dual bounds above are taral's own numbers (not recomputed by independent code); the rows stay FAIL under the stated rule.
No solver change was needed for these five. Not run here: the 18 rows without a HiGHS reference and the 4 numerical failures.
