# R4 final classification - MIPLIB3 65/65 at c4800e9 (Oct 8, 2026)

Protocol: v3 runner on Azure D4s_v4 (co-located, 2 workers, isolated dirs), taral --time-limit 300
--json --sol; HiGHS 1.15.1, 1 thread, 300 s reference. misc06's reference rechecked at
mip_rel_gap=0 after the false-WRONG resolution (engine independently confirmed; raw TARAL point
max violation 6.39e-13); the other 64 references are as-run HiGHS default-gap. All 65 original models hash-matched; ALL reported points pass the 1e-6
original-model audit; all available optimal-reference bound checks pass (engine point-audit JSON).
Audit scope note: 61 of 65 cases returned points; 4 have check=null (no incumbent) - the audit
statement is 'all REPORTED points pass', never 'all 65 raw points'.

Final ledger (65/65, corrected):
- 31 proven_agree: taral proved optimal (gap 0) and matches the HiGHS reference within 1e-6
  relative (misc06's reference rechecked at zero gap; the rest as-run default gap, agreement
  within 1e-6 in all cases). misc06 included (published-optimum match, reference artifact resolved).
- 5 value_reached_unproven: 10teams, air04, nw04, qiu, vpm1 - incumbent equals the reference
  optimum within 1e-6 (worst 9.1e-16) but the 300 s cap hit before the proof closed. "Unproven",
  never agreement.
- 17 differs (neutral label, was v2 "incumbent_below_ref"): unproven mid-search incumbents at the
  cap differing from the reference optimum: cap6000, fiber, fixnet6, gesa2, gesa2_o, gesa3, harp2,
  mod011, modglob, noswot, p0548, p2756, pp08a, pp08aCUTS, rout, set1ch, vpm2. None claims optimal.
- 12 other: no taral incumbent at cap (gesa3_o, mitre) or no proven reference (both solvers capped:
  arki001, dano3mip, danoint, fast0507, markshare1, markshare2, mas74, mkc, seymour, swath).
  Neither agreement nor disagreement; labelled, never dropped.
- 0 WRONG_DISAGREE.

Value-reached total: 36/65 (31 proven + 5 unproven). Oct 5 comparison: 28/65 proven, 35/65
value-reached, mixed builds; this run is a single clean c4800e9 build, one protocol, so 31/65
proven / 36/65 value-reached is the clean comparable figure.

Known-defect watch: zero optimal-claiming rows disagree with the reference. The near-integrality
defect (found Oct 8 via the R2 corpus pattern, confirmed in mainline by the engine) has no field
instance in this corpus; results still PREDATE its fix (old-baseline run by design). Engine's
standing caveat: not an exact global proof beyond the numeric gate while the near-integrality fix
is separate.

air04: byte difference vs the engine's older cache is ONLY the MPS comment line (BEST SOLN 56137
vs 56138); parsed model field-identical. Caveat adopted: comment-published "best" values are never
used as verified optimum references.

Runner v4: mip_rel_gap=0 pinned on the reference (fixed in my master m3run_vm.py).

## Suggested README reword (R4)
Current: - [ ] **R4 - Mixed-integer search:** prototype [branch and bound with propagation](src/milp.cpp); full cutting-plane, presolve and heuristic coverage planned.

Proposed (assuming results land at results/r4_miplib3_20261008/):
- [x] **R4 - Mixed-integer search:** prototype [branch and bound with propagation](src/milp.cpp); [measured on all 65 MIPLIB 3 instances](results/r4_miplib3_20261008/) (300 s cap, HiGHS 1.15.1 reference with misc06 rechecked at zero gap, every reported point independently feasibility-audited at 1e-6): 31/65 proven optimal matching the reference, 5 more reached the reference optimum value without completing the proof, 29 unproven at the cap, 0 wrong answers on the tested corpus; a near-integrality correctness defect found after this run is being repaired separately, and these results predate that fix; full cutting-plane, presolve and heuristic coverage planned.

Honesty notes on the line: the 29 unproven (17+12) are in the line itself; "0 wrong answers" is
scoped to the tested corpus at the 1e-6 gate; the defect caveat is in the sentence, not buried.

## Repo-lane wording corrections adopted (Oct 8 03:25)
1. Audit says "all reported points pass" (61/65 returned points; 4 have check=null/no incumbent).
2. Reference is as-run HiGHS default gap; ONLY misc06 was rechecked at mip_rel_gap=0. Never call
   the whole reference exact-gap. Repo preserves the original 30+1-false-WRONG ledger plus a
   clearly labelled corrected ledger.
