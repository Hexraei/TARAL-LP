# Post-Q-merge second-host Netlib confirmation

Kaggle replay of source `a240203` with the integrated Q patches: 93/93 passing
in each mode when pilots are counted, strict 92/93. Only `greenbea` is nonstrict;
no case fails. The CSVs are independently compared with the `58f77af` evidence:
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
