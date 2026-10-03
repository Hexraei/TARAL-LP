# TARAL-LP re-measure: Kennington / Mittelmann LP benchmark instances, merged engine

**Not the Netlib 93-gate.** Different instance set, 300 s solver cap / 600 s HiGHS cap, different host. None of these rows enter the 93 denominator.
Engine: origin/main 58f77afbb8e85e51ad2cb9971be9e12d8e2bd96b, g++ 13.3.0 -O2 (zero diagnostics under -Wall -Wextra -Wpedantic). HiGHS 1.15.1 defaults, 600 s cap.
Host: Intel Xeon @ 2.80 GHz VM, 4 cores, 16 GB, Linux 6.18.44. The Oct 2 baseline ran on a previous CPU host; wall times are not comparable across hosts. Details: HOST.txt.
Each case was run once per mode for the table; the six rows whose status changed vs Oct 2 were run a second time. Stages ran one at a time. Every optimal row was checked by the independent checker (benchmarks/orig_check.py, 1e-6 relative).

Summary of changes vs Oct 2 (33 cases, 300 s cap):
1. Wrong answers: 0. All 40 optimal results (21 dual + 19 default) are checker-feasible and match HiGHS within 1e-6 relative (worst objective error 7.9e-15, worst row violation 2.7e-12, worst bound violation 6.5e-11).
2. Dual (--method dual): 21 optimal / 12 time_limit, versus 19 / 14 on Oct 2. Newly optimal: fome12 (273.3 s, repeat 261.5 s) and fome21 (250.2 s, repeat 251.9 s). Newly timed out: none. Mismatches: none.
3. Default method: 19 optimal / 14 time_limit, versus 15 / 18 on Oct 2. Newly optimal: fome11 (137.7 s, repeat 136.4 s), cre-b (86.5 / 87.3 s), cre-d (79.4 / 78.3 s), pds-20 (122.5 / 121.0 s). Newly timed out: none. Mismatches: none.
4. The default-method gains come from the merged primal-budget fallback: after about 60 s of stalled primal simplex the engine hands the rest to the dual simplex (recorded in each row's message). 7 of the 19 default optimal rows used that route (fome11, rail582, ken-13, cre-b, cre-d, pds-10, pds-20); 12 finished in primal simplex alone. The default method is therefore not a pure-primal measurement any more.
5. The dual-mode gains are not shown to be an engine effect: dual iteration counts match Oct 2 where both finished (e.g. ken-07 2929, ken-13 31585), and wall times on this host are about 1.2x to 2.0x shorter than on the previous CPU host for dual rows above 5 s (fome11 152.6 s to 74.7 s). fome12 needs 86.8k dual iterations; Oct 2 reached 57.9k at its cap. The pre-merge engine was not run on this host, so host and engine cannot be separated.
6. Current implementation times out (300 s cap) on 12 cases in dual mode: fome13, rail2586, rail4284, ken-18, pds-30..100. Default method times out on those 12 plus fome12 and fome21.
7. HiGHS finished 31 of 33 within 600 s; rail2586 and rail4284 hit its cap, same as Oct 2. So 2 time_limit rows have no reference, and neither engine mode produced an answer there.
8. HiGHS is faster than the engine on every row both finished: 21 of 21 dual-optimal rows and 19 of 19 default-optimal rows. No speed parity or advantage is claimed.
9. Repeats: all 6 rows reproduced status, objective and checker values; wall times within about 5%; iteration counts are identical in dual mode and differ by up to about 1.5% in default mode.
10. Provenance: the expanded MPS sha256 matches the Oct 2 table for 33 of 33 cases; download sha256 for 32 of 33 (pds-20: I fetched netlib's copy; plato's copy has the Oct 2 hash, expanded MPS identical).

Notes: all inputs are emps-encoded ("compressed MPS") and were expanded with netlib's emps.c. netlib's kennington directory holds only pds-02/06/10/20, so pds-30..100 come from plato lptestset/pds (URLs per case in the CSV). A time_limit row is "no answer", not a correctness failure. Wall is process wall including parse; HiGHS wall is solver run() only (read time is a separate column), which favours HiGHS in the comparison. Rows/cols/nnz are copied from the Oct 2 table (identical MPS bytes). The second optional batch was not started.
