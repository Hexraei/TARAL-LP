# Post-Q-merge second-host Netlib confirmation

## Reading the technical terms

LP means linear programming; MILP means mixed-integer linear programming; QP means quadratic programming. IPM is an interior-point method. PDHG is primal-dual hybrid gradient, the approximate GPU method. KKT (Karush-Kuhn-Tucker) checks test feasibility and optimality conditions. A "gate" is the stated validation rule, not an exact-arithmetic proof. "Strict" means the separately stated tighter tolerance; a non-strict pass meets the ordinary rule but not that tighter check. A ledger is a per-case result table. Source hashes, file paths and command flags are retained only so engineers can reproduce a measurement. Historical measurements are not current-source claims.

JSON is the machine-readable result format. Warm starts reuse a previous solution or simplex basis; cold starts do not. fp64 is double-precision floating point; FMA is fused multiply-add.


Kaggle replay of evaluated source snapshot `a240203` with the integrated Q patches: 93/93 passing
in each mode when pilots are counted, strict 92/93. Only `greenbea` is nonstrict;
no case fails. The CSVs are independently compared with the evaluated source snapshot `58f77af` evidence:
zero changes in case identity, verdict, passed, strict, engine status or printed
objective in both modes. Wall times differ and are not a speedup claim.

`summary.txt` and the logs retain the harness's original pilot-excluded display
(91/93, strict 90/91). The raw CSV rows show both pilots passing; including them
produces 93/93 and strict 92/93. No raw ledger or log was changed to make this
count. The summary's wording about a "headline convention" is historical;
the current README headline uses pilots counted.

Notebook identifier, archive/source fingerprints, flags, protocol and artifact
hashes are in `provenance.json`. The source archive was verified by the executor;
this record does not contain a fresh independent full QP harness run. Equality
harness relabel code is not exercised by this Netlib run. Printed objectives
match byte-for-byte as CSV strings; original in-memory double bytes are absent.
