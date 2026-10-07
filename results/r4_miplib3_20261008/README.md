# R4: 65-instance MIPLIB 3 old-baseline measurement

Measured engine: c4800e9944587f221cfe25b9aebfa2a31750ff8b. Azure D4s_v4,
two colocated workers, 300 s per solver case. This is not a speed-comparison claim.
The run predates the separate near-integrality correctness fix. No claim of an
exact global proof beyond the 1e-6 numerical gate is made.

## Corrected classification

- 31 solver-optimal at the numerical gate, objectives match the reference within 1e-6 relative.
- 5 additional value matches without closing the proof: 10teams, air04, nw04, qiu, vpm1.
- 17 differing unproven incumbents; 12 other cases without a TARAL point or a proven reference.
- 0 wrong optimal claims on this tested corpus after the misc06 reference correction.

There are 65 audited cases but only 61 reported points. Four cases have no point:
arki001, gesa3_o, mitre, swath. All 61 returned points pass the independent
original-model feasibility audit at 1e-6. All 55 available reference-bound checks
pass; the other ten have no proven reference.

"Optimal" here uses TARAL's numerical gap tolerance, not gap exactly zero:
bell3a, bell5, dcmulti, egout, mas76 and pk1 have nonzero gaps below 1e-6.
The supplied final-classification memo's phrase "gap 0" is therefore not accurate
for every solver-optimal row. Raw outputs are preserved to make this visible.

## Reference correction and provenance

The as-run reference uses HiGHS 1.15.1 defaults, 1 thread, 300 s. It is NOT an
all-case zero-gap reference run. Only misc06 was rerun with mip_rel_gap=0:
its objective 12850.860737382538 is identical to TARAL's, whereas the default-gap
reference stopped at 12851.076291564139 (gap 7.36e-05). The tight rerun has bound
equal to objective and zero gap. Its model hash matches the ledger and raw model.
`ledger_corrected.csv` applies this one correction and renames the 17 neutral
`incumbent_below_ref` labels to `incumbent_differs_from_ref`. No solver output is changed.
The correction's wall time was not supplied; that field is blank, not borrowed
from the original run.

`as_run/` preserves every supplied original ledger, reference, output and log,
including the original misc06 WRONG_DISAGREE and original summary. `point-audit.json`
is the independent engine audit as supplied (its misc06 reference objective is
still the as-run value). `misc06-tight-reference.json` supplies the correction.
`provenance/` contains runners v3/v4, build/source manifests and the instance list.
V4 pins zero reference gap for future runs, not this past 65-instance run.

All 65 `mps/` files were independently hash-checked against the original ledger.
air04's comment-line BEST SOLN is not an optimum reference. Only solver references
are used. The supplied provenance memo's archive inventory overstates its contents:
the point-audit and tight-reference JSON arrived separately, and are stored here
separately, not represented as inside the original archive.

Original results archive SHA256: a1d647958a85d7314e90e3fdedea781a6eab5d04f8bedd785eef58c60190e6cd
Original models archive SHA256: f76b320f3b62eb783c132ad9bda4e5a1bed8ec751035eea1e4b48366bb2bf586
No new solver benchmark was run while packaging this evidence.
