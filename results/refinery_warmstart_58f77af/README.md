# Synthetic refinery warm-start measurement

39 what-if LP rows across three synthetic refinery-shaped bases. All cold and warm checker labels are `ok` at the reported 1e-6 oracle tolerance; the recorded wrong count is zero. Independently recomputed median warm/cold iteration ratio: **0.3022847100**. This is an iteration-count result, not a wall-time speedup claim.

Engine: `58f77af`, default simplex. Scripts: main `6c621e0`, [benchmark driver](../../benchmarks/refinery_warmstart.py). Oracle: SciPy 1.17.1 `linprog(method="highs")`, bundled HiGHS 1.12.0. Host and binary metadata are recorded in [provenance](provenance.json).

`warm_used` is inferred from absence of fallback text in stderr, not positive evidence of the warm-start route. Wall times are 3-19ms and are too short for a stable speed comparison. These are synthetic LPs, not Netlib cases or operational plant data. The CSV carries checker labels, not complete point/certificate artifacts.

[Raw table](warmstart_table.csv) · [Raw summary](warmstart_summary.json)
