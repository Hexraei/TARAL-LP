# LP certificate quality and solve acceptance

Strict original-model KKT quality is reported separately from solve acceptance.
`certificate_quality` is `pass`, `fail`, or `unknown`, using the unchanged 1e-8
primal/dual/gap/complementarity metric. It is not an acceptance veto. Unknown is
reported for unavailable/nonfinite evidence or an unresolved/nonoptimal solve.
The baseline internal 1e-7 original-model feasibility check remains unchanged.
The benchmark's original-model 1e-6 row/bound rule remains unchanged too.

GREENBEA can therefore be `optimal` with quality `fail`; it does not meet the
strict 1e-8 metric. This is explicit, not hidden or described as a strict pass.
The independent checker still exits nonzero for strict fail, independently of
whether the solve meets canonical acceptance.

Existing `primal_res`, `max_row_viol` and other labels retain their definitions.
New original-coordinate per-row arrays report:

- `row_violation_abs`: max(0, lower - activity, activity - upper).
- `row_term_magnitude`: sum_j |a_ij x_j| + |rhs_i|. For equality rows RHS is the
  unique bound. For violated one-sided/ranged rows RHS is the violated finite
  endpoint. For feasible non-equality rows RHS is zero, unless the activity is
  within 1e-12 * (sum_j |a_ij x_j| + |endpoint|) of a finite endpoint; then that
  endpoint is the RHS, whether the activity lands just inside or just outside it.
  Reason: summation order or fused multiply-add changes the last bit of the
  activity, and a rule that switches on the sign of a 1e-15 violation made the
  magnitude jump by |endpoint| (random_13: engine -7.500000000000001 vs checker
  -7.5 against a lower bound of -7.5, magnitudes 37.46 vs 29.96). The band is
  `kEndpointBand` in `src/simplex.cpp` and `ENDPOINT_BAND` in
  `benchmarks/simplex_certificate_check.py`; keep them equal. Raw violations,
  violated endpoints and equality rows are unchanged, and the checker still
  flags any real mismatch beyond its 1e-8 relative tolerance.
- `row_violation_magnitude_scaled`: raw absolute violation / row term magnitude.
  If the denominator is zero, scaled violation equals the raw violation.
- `max_row_violation_magnitude_scaled`: maximum of the per-row scaled array.

These are reporting-only. They never replace the canonical RHS-normalized rule.
A tiny violated zero-RHS single-term row can have scaled violation 1 even with
very small raw violation. Large cancellation rows can have tiny magnitude-scaled
violation while still failing strict RHS-normalized quality. Both are honest.

## Acceptance policy versus earlier commits

Canonical acceptance (the fixed rule in the README) is compared against two
earlier states. The two comparisons differ and must not be merged into one claim.

- Versus `90430a4` (the documented baseline): canonical acceptance is unchanged.
  GREENBEA was already a pass there, with the strict check failing.
- Versus `eac35f2` (an earlier version of this branch that vetoed a solve whose
  strict quality was `fail`): acceptance CHANGED. The stricter quality veto was
  removed on purpose. GREENBEA is now accepted as `optimal` with the same
  objective (-72555248.1298), `certificate_quality` = `fail` and `primal_res`
  about 1.44e-8 against the strict 1e-8. The pass counts rise by one in each
  method (primal 89 to 90, dual 90 to 91 of 93 cases). The strict counts
  (89/90 primal, 90/91 dual) are unchanged, so the strict check still fails
  GREENBEA and this is a non-strict pass, not a strict one.

Evidence for the `eac35f2` figures is the independent reviewer's report; that
commit is not in this checkout, so those counts were not re-measured here.
Re-measured here at the current code with the pinned GREENBEA file (sha256
`3d978068...fdf`): primal `optimal`, objective -72555248.12984599, 9037 iterations,
quality `fail`, `primal_res` 1.441e-8; dual `optimal`, 5640 iterations, quality
`fail`, `primal_res` 6.14e-8. `certificate_quality_tests.py --greenbea` passes.

```
g++ -O3 -std=c++17 -Wall -Wextra src/*.cpp -o /tmp/taral
python3 benchmarks/simplex_certificate_tests.py --engine /tmp/taral
g++ -O2 -std=c++17 -Wall -Wextra benchmarks/certificate_quality_api_tests.cpp src/{mps,simplex,dual,lu}.cpp -o /tmp/quality-api
/tmp/quality-api
```
