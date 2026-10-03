# Kennington and Mittelmann: 33-case measurement record

Measured evidence with a scoped same-host comparison. This is a different benchmark selection from the 93-case
Netlib validation, with a 300-second engine limit and 600-second HiGHS limit.
The evaluated source snapshot is `58f77af`; it predates the later quadratic
programming changes. `ledger.csv`, `SUMMARY.md`, `HOST.txt` and
`inputs_manifest.csv` are retained in their corrected original CRLF byte form. The first archived copies used LF line endings; the field values are unchanged. Artifact hashes now match the corrected delivery. `HOST.txt` remains unchanged. `SUMMARY.md` items 5 and 9 were later revised from the full measurement report; its original and revised hashes are recorded in provenance.

Independent table checks find 21 optimal / 12 time-limit results using dual
simplex and 19 / 14 using the default method. All recorded optimal results
match the reference under the stated checks. Seven of the 19 default optimal
results switched from primal to dual simplex; the default column is not a
pure-primal measurement. `primal_route` records this method choice.

The same-host comparison in `isolation_442ca16_vs_58f77af_fome12_fome21.csv`
supports an engine-version effect on two dual-simplex cases: the earlier evaluated
source snapshot `442ca16` reaches the 300-second limit on both, while `58f77af`
finishes `fome12` in 270.56 seconds and `fome21` in 253.76 seconds. The checked
relative objective errors are 7.9e-15 and 0. Iteration throughput is about 1.25x
higher for the newer build. These are single runs per build/case on the same
reported host, not a repeated performance study or an isolated attribution to
one individual patch. The attached comparison CSV is a transcription of the
measurement report, not the executor's original file bytes; its own hash is
recorded separately.

This supersedes the unresolved attribution in the original `SUMMARY.md` item 5 and adds
the single-run isolation context to item 9 for these two cases only. The four
newly optimal default-method cases still have no same-host control. Seven of
the 19 default optimal results switched to dual simplex, so their improvement
must not be described as pure-primal progress. Cross-host wall times for the
full selection remain confounded by the 1.2-2.0x faster host.

HiGHS is faster on every case where both finish. Neither solver finishes
`rail2586` or `rail4284` within its own limit, so those have no reference answer.

The raw summary records six changed statuses and their repeats. All are retained
as measurements, not a claim that future runs or every industrial model will
behave the same. Manifest MPS hashes independently match the result table for
33/33 cases; original instance bytes and full per-case logs are not included
in this delivery. Artifact hashes and checked counts are in `provenance.json`.
