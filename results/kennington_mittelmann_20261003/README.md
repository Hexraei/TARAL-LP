# Kennington and Mittelmann: 33-case measurement record

Raw evidence only. This is a different benchmark selection from the 93-case
Netlib validation, with a 300-second engine limit and 600-second HiGHS limit.
The evaluated source snapshot is `58f77af`; it predates the later quadratic
programming changes. `ledger.csv`, `SUMMARY.md`, `HOST.txt` and
`inputs_manifest.csv` are retained byte-for-byte as received.

Independent table checks find 21 optimal / 12 time-limit results using dual
simplex and 19 / 14 using the default method. All recorded optimal results
match the reference under the stated checks. Seven of the 19 default optimal
results switched from primal to dual simplex; the default column is not a
pure-primal measurement. `primal_route` records this method choice.

The apparent dual improvement versus the previous host is not yet attributable
to solver changes. The reported host is 1.2-2.0x faster; a same-host comparison
with the earlier source is pending. No improvement headline is published here.
HiGHS is faster on every case where both finish. Neither solver finishes
`rail2586` or `rail4284` within its own limit, so those have no reference answer.

The raw summary records six changed statuses and their repeats. All are retained
as measurements, not a claim that future runs or every industrial model will
behave the same. Manifest MPS hashes independently match the result table for
33/33 cases; original instance bytes and full per-case logs are not included
in this delivery. Artifact hashes and checked counts are in `provenance.json`.
