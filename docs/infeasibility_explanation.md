# Verified LP infeasibility explanation

`taral MODEL.mps --explain-infeasible OUT.json` explains why the LP relaxation of a model is
infeasible and reports a smallest-weight relaxation that would make it feasible. Integrality
and the objective are ignored (the solve uses zero cost). Models with a quadratic objective are
not special-cased: the objective is not used either way.

## Definitions

- **Irreducible (this tool's meaning):** row-irreducible under retained column bounds. The
  explained row set, together with all column bounds, has a verified Farkas certificate, and
  for every row in the set, removing that row leaves a system with a stored witness point that
  satisfies all remaining rows and all column bounds at relative tolerance 1e-7.
  It is not a globally minimum-cardinality infeasible subsystem. Different irreducible sets can
  exist, and the one reported depends on the deletion order (least Farkas multiplier first).
  Column bounds are never dropped, so a bound that takes part in the conflict is listed under
  `bound_columns` (with its multiplier) but is not a removable element.
- **Relaxation report:** the optimum of an elastic LP over all rows. Each row gets non-negative
  slacks on its lower and upper side, weight `1/(1+max finite |row bound|)`, column bounds
  retained. It minimises the weighted L1 sum of the slacks. The optimum value is unique; the
  point and the set of relaxed rows need not be. It does not minimise the number of relaxed rows.
- **Statuses:** `irreducible`, `reduced_unproven` (some row's removal test was inconclusive, listed
  in `unproven_rows`), `bounds_only` (a column has lower > upper), `relaxation_feasible` (the LP
  relaxation is feasible, nothing to explain), `no_verified_proof` (no certificate passed
  verification; exit code 5).
- Certificates are floating-point, verified at the engine's existing margin (`docs/certificates`
  semantics, 1e-8); they are not exact-rational proofs.

## Behaviour notes (after independent review)

- The time limit covers every phase. The elastic solve is skipped (`relaxation.status` =
  `not_run_time_limit`) when no budget is left after the deletion phase; no minimum one-second grace.
  The row proof status is then reported without a relaxation report. Timed-out removal tests leave
  rows unproven (`reduced_unproven`).
- An unwritable output path exits 5 with a message on stderr and stdout; no success is claimed.
- The weighted L1 side relaxation is reporting-only and is NOT claimed minimal. Its elastic costs are multiplied by one exact
  power of two `scale_factor` s = 2^k, k = max(ceil(log2(1/smallest positive weight)), 20), capped at 30 (a non-finite
  log2(1/weight), from an overflowing 1/weight, clamps to 30). Reason: the smallest weights (about 3e-8 on Netlib cplex1) are within a
  factor of about 30 of the simplex's absolute 1e-9 reduced-cost tolerance, and the solve can stop early at a worse point. Relative weights
  and the optimal set are unchanged and the division by s is exact. `relaxation.objective` is in original units (solver value / s).
  `scale_factor` is null only when the scaled LP was never attempted.
- `relaxation.certificate_quality` is the gated solver's own strict 1e-8 KKT label of the SCALED LP (`certificate_quality_model` =
  `scaled_lp`); it is not an exact proof. Values: pass | fail | unknown from the solver; `unverified_point` keeps the solver label and has a
  null objective; `not_available` when the scaled LP was attempted but produced no label (non-optimal solve); `not_run` when it was not
  attempted. `kkt_gap` and `kkt_gap_scaled_model` are the scaled model's relative gap; `kkt_gap_abs_unscaled` = |primal - dual objective| / s
  in original units. A "fail" can come from the absolute 1e-8 primal residual alone (Netlib cplex1, native build: 1.49e-8) and is
  unchanged by scaling; the gate is not loosened.
- `relaxation.max_violation` (original units, long double recheck, all below the 1e-7 witness tolerance). Figures and who measured them,
  builds are not blended: cplex1 native unpatched 5.05e-9 comes from the author's elastic bundle on 692dac0d (native -O3 -march=native build,
  explain.json SHA256 d53528a27ffbae4d24f6cde9e014a764484a495a014a9bbf01026762a810997f), which the verification side independently verified.
  The author's own runs on 692dac0d (generic -O2 build for cplex1, native build for qual; before = unpatched main, after = this patch):
  cplex1 generic 1.609e-8 -> 2.290e-8; qual native 4.75e-12 -> 5.44e-12. Whether the verification side measured these same figures on
  9b09e63 is not established here (matching numbers do not establish an identical tree), so no 9b09e63 attribution is made for them.
- Measured on Netlib cplex1 and qual (netlib.org/lp/infeas, EMPS-expanded; local runs, generic -O2 -march=x86-64 -mtune=generic
  -fno-fast-math -ffp-contract=off and native -O3 -march=native builds): cplex1 generic objective 0.38937931 -> 0.37766446 (s = 2^25),
  qual native 0.00541015975 -> 0.00541015959 (s = 2^20). The bracketing analysis (scripts and results, diagnosis pack SHA256
  44a4288789b02f5b72e70e46b6594c23cbd1d39badc0041d3ccc4afdd15f07b9) was produced by the independent review side, not by this repository:
  a conditional weak-duality lower bound and an exact-rational feasible upper bound bracket cplex1 .377664455209 / .377664551 and qual
  .005410159586 / .005410159615; the bounds are conditional on artificial boxes on infinite-bound columns. A true optimum is not proved here.
- The checker's relaxation oracle (HiGHS through scipy) multiplies its objective by the fixed power of two 2^20 and divides back: its
  absolute dual tolerance otherwise stops its simplex early on small weights (qual: 0.00542051 instead of 0.00541016). The checker also
  recomputes the objective of the reported point exactly (Fractions from the parsed binary doubles, x not clipped, exact column-bound
  violation reported separately) and compares it to the reported value within 1e-6 * (1 + |objective|); the scaled float oracle is compared
  within the same 1e-6 * (1 + |objective|) and is an independent second check.
- Duplicate coefficients in an MPS COLUMNS section are summed by the engine. HiGHS keeps the first one, so the
  checker normalizes the model (sums duplicates, free format only) before using HiGHS as the oracle.
- The checker never trusts `certificate_verified`: it re-checks multiplier dimensions, finiteness, signs,
  stationarity (<= 1e-9 scaled by the L1 norm) and a contradiction margin > 1e-8 on the original model,
  rejects NaN/non-finite JSON, wrong-length or non-finite witnesses and bad indices, requires oracle
  status infeasible (not merely non-optimal) and fails if the oracle read or solve fails.

## Measured (this branch)

Measured by `benchmarks/infeasibility_explanation_tests.py` and the independent checker
`benchmarks/infeasibility_explanation_check.py`:

- `benchmarks/refinery_stress/infeasible_supply_lp.mps`: 26 explained rows, status `irreducible`,
  31 LP solves, under 0.01 s on the measuring host. The checker (rows parsed by HiGHS' reader,
  feasibility via the HiGHS LP solver as test oracle) confirms the 26 rows are infeasible, that
  dropping each row alone is feasible, that every witness is feasible at 1e-7, and that the
  relaxation point and objective match the oracle optimum within 1e-6 * (1 + |objective|) (the oracle objective scaled by 2^20 and divided back,
  and the reported objective recomputed exactly from the point, see above).
- Small hand-built cases: two-row conflict, bounds conflict, feasible model, and three corrupted
  explanations (dropped row, bad witness, wrong objective) that the checker rejects.
- The feasible refinery fixtures report `relaxation_feasible`.

## Inherited / attributed (not measured here)

- The full Farkas support of the same fixture was 54 of 72 rows in an earlier run of the
  pre-existing certificate path; HiGHS 1.15.1 `iis_strategy=2` returned 37 rows. Neither count
  is evidence of irreducibility or minimality (the HiGHS documentation does not claim it), and
  the tool's 26 rows are not shown to be smaller than a minimum-cardinality set.
- Behaviour of other solvers' IIS tools is attributed to their documentation only.

## Not run / not claimed

- No dedicated infeasible corpus sweep (MIPLIB/Netlib); only 150 generated random LPs (seeds 7000-7149); the historical 26-time-limit MIPLIB ledger is not a
  baseline for this feature.
- Irreducibility is not shown for integer models (LP relaxation only).
- No exact-arithmetic or rational certificate; no minimum-cardinality or minimum-row-count search.
- No timing comparison against HiGHS IIS.
- Budget behaviour is tested with a simulated clock after the deletion phase, not with a real overrun.
- Fixed-format MPS with spaces in names is not supported by the checker's duplicate normalizer.
