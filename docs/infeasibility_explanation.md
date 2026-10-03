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

## Measured (this branch)

Measured by `benchmarks/infeasibility_explanation_tests.py` and the independent checker
`benchmarks/infeasibility_explanation_check.py`:

- `benchmarks/refinery_stress/infeasible_supply_lp.mps`: 26 explained rows, status `irreducible`,
  31 LP solves, under 0.01 s on the measuring host. The checker (rows parsed by HiGHS' reader,
  feasibility via the HiGHS LP solver as test oracle) confirms the 26 rows are infeasible, that
  dropping each row alone is feasible, that every witness is feasible at 1e-7, and that the
  relaxation point and objective match the oracle optimum (relative 1e-6).
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

- No MIPLIB or Netlib infeasible-set sweep; the historical 26-time-limit MIPLIB ledger is not a
  baseline for this feature.
- Irreducibility is not shown for integer models (LP relaxation only).
- No exact-arithmetic or rational certificate; no minimum-cardinality or minimum-row-count search.
- No timing comparison against HiGHS IIS.
